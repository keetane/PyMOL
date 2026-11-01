"""PyMOL commands for retrieving protein structures.

Load this file in PyMOL with::

    run chemistry.py

Then retrieve a human AlphaFold model by gene/protein name or UniProt accession::

    af2 TP53
    af2 P04637
    af2 "cellular tumor antigen p53"

Print the SMILES for HETATM records in ``sele`` (or another selection) and
copy it to the system clipboard::

    printsmiles
    printsmiles ligand

Build and load a 3D object from a SMILES string::

    smiles "CCO"
    smiles "CCO", ethanol

Search PubChem by name, copy its SMILES, and load its 3D SDF conformer::

    pubchem aspirin
    pubchem "acetylsalicylic acid"

Build a 3D peptide from a one-letter sequence and an optional cyclization mode::

    peptide ACDE
    peptide ACDE, 1, cyclic_ACDE

Create a directory relative to PyMOL's current working directory::

    mkdir results/models

Delete one file or files matching a non-recursive glob pattern::

    rm obsolete.pdb
    rm *.cif

Build an aligned structure session for a human protein target::

    target JAK1

Print docking scores for the currently displayed state of enabled objects::

    printscore
    printscore docked, CNNscore

The downloaded PDB file is saved in PyMOL's current working directory and is
loaded as a new object.  Only Homo sapiens (NCBI taxonomy ID 9606) entries are
accepted.
"""

from __future__ import annotations

import csv
import datetime as dt
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, urlparse
from urllib.request import Request, urlopen


def _reinvent_safe_filename(value: Any) -> str:
    result = re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_")
    return result[:64] or "UNKNOWN"


def _reinvent_detect_device() -> str:
    import torch

    if torch.cuda.is_available():
        try:
            tensor = torch.ones((4, 4), device="cuda")
            torch.cuda.synchronize()
            if (tensor @ tensor).device.type == "cuda":
                return "cuda"
        except (RuntimeError, AssertionError) as exc:
            print(f"[REINVENT] CUDA test failed: {exc}", flush=True)

    backend = getattr(torch.backends, "mps", None)
    if backend is not None and backend.is_available():
        try:
            tensor = torch.ones((4, 4), device="mps")
            result = tensor @ tensor
            torch.mps.synchronize()
            if result.device.type == "mps":
                return "mps"
        except (RuntimeError, AssertionError, TypeError) as exc:
            print(f"[REINVENT] MPS test failed: {exc}", flush=True)
    return "cpu"


def _reinvent_reported_device(log_text: str) -> str:
    matches = re.findall(
        r"Using (CPU|CUDA|MPS)\b[^\n]*",
        log_text,
        flags=re.IGNORECASE,
    )
    return matches[-1].strip() if matches else "not reported"


def _reinvent_worker(manifest_file: str) -> None:
    from rdkit import Chem

    manifest_path = Path(manifest_file).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    job_directory = manifest_path.parent
    manifest["device"] = _reinvent_detect_device()
    manifest["device_report"] = {}
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[REINVENT] Auto-selected device: {manifest['device']}", flush=True)

    raw_molecule = Chem.MolFromMolFile(
        str(job_directory / "selected_ligand.mol"),
        removeHs=False,
        sanitize=True,
    )
    if raw_molecule is None:
        raise RuntimeError(
            "Ligand could not be parsed as a chemically valid MOL. "
            "Verify bond orders and charges in PyMOL."
        )
    molecule = Chem.RemoveHs(raw_molecule)
    input_smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
    if not input_smiles or Chem.MolFromSmiles(input_smiles) is None:
        raise RuntimeError("Could not obtain a valid SMILES from selected ligand.")
    (job_directory / "input.smi").write_text(input_smiles + "\n", encoding="utf-8")
    (job_directory / "input_smiles.txt").write_text(
        input_smiles + "\n",
        encoding="utf-8",
    )

    executable = Path(manifest["runtime"]["executable"])
    if not executable.is_file():
        raise RuntimeError(f"Missing REINVENT executable: {executable}")

    merged_csv = Path(manifest["output_csv"])
    partial_csv = merged_csv.with_suffix(".partial.csv")
    header = None
    total_rows = 0
    for prior_value in manifest["priors"]:
        prior = Path(prior_value)
        tag = _reinvent_safe_filename(prior.stem)
        target_csv = job_directory / f"{tag}.csv"
        config_path = job_directory / f"{tag}.toml"
        config_path.write_text(
            'run_type = "sampling"\n'
            f'device = {json.dumps(manifest["device"])}\n\n'
            "[parameters]\n"
            f'model_file = {json.dumps(str(prior))}\n'
            f'smiles_file = {json.dumps(str(job_directory / "input.smi"))}\n'
            f'output_file = {json.dumps(str(target_csv))}\n'
            f'num_smiles = {manifest["num_smiles"]}\n'
            f'temperature = {manifest["temperature"]}\n'
            'sample_strategy = "multinomial"\n'
            "unique_molecules = true\n"
            "randomize_smiles = true\n",
            encoding="utf-8",
        )
        print(
            f"[REINVENT] Starting {prior.name} on {manifest['device']}",
            flush=True,
        )
        run_environment = os.environ.copy()
        if manifest["device"] == "mps":
            run_environment.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
        result = subprocess.run(
            [
                str(executable),
                "--device",
                manifest["device"],
                str(config_path),
            ],
            cwd=str(job_directory),
            env=run_environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            check=False,
        )
        output_text = result.stdout or ""
        (job_directory / f"{tag}.log").write_text(output_text, encoding="utf-8")
        reported = _reinvent_reported_device(output_text)
        manifest["device_report"][tag] = {
            "requested": manifest["device"],
            "reinvent_log": reported,
            "exit_code": result.returncode,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        print(
            f"[REINVENT] Device: requested={manifest['device']}, "
            f"REINVENT log={reported}",
            flush=True,
        )
        if manifest["device"] == "mps" and reported.lower().startswith("using cpu"):
            print(
                "[REINVENT] WARNING: REINVENT reports CPU despite explicit "
                "--device mps. GPU computation is not confirmed.",
                flush=True,
            )
        if result.returncode != 0 or not target_csv.is_file():
            print(
                f"[REINVENT] FAILED {tag}: exit={result.returncode}. "
                f"See {tag}.log",
                flush=True,
            )
            continue

        with target_csv.open(newline="", encoding="utf-8-sig") as input_handle:
            reader = csv.DictReader(input_handle)
            fields = reader.fieldnames or []
            if not fields:
                print(f"[REINVENT] Empty CSV for {tag}", flush=True)
                continue
            if header is None:
                header = [field for field in fields if field != "prior"] + ["prior"]
                with partial_csv.open("w", newline="", encoding="utf-8") as output_handle:
                    csv.writer(output_handle).writerow(header)
            count = 0
            with partial_csv.open("a", newline="", encoding="utf-8") as output_handle:
                writer = csv.DictWriter(
                    output_handle,
                    fieldnames=header,
                    extrasaction="ignore",
                )
                for item in reader:
                    item["prior"] = prior.stem
                    writer.writerow(item)
                    count += 1
            total_rows += count
        print(f"[REINVENT] Finished {tag}: {count} rows", flush=True)
        if partial_csv.exists():
            shutil.copyfile(partial_csv, merged_csv)

    if header is None:
        raise RuntimeError(
            f"All model runs failed. Inspect log files under {job_directory}"
        )
    print(f"[REINVENT] Complete: {merged_csv} ({total_rows} rows)", flush=True)


if len(sys.argv) == 3 and sys.argv[1] == "--reinvent-worker":
    _reinvent_worker(sys.argv[2])
    raise SystemExit(0)


from pymol import cmd, preset, util


_AFDB_BASE = "https://alphafold.ebi.ac.uk"
_AFDB_SEARCH_URL = f"{_AFDB_BASE}/api/search"
_AFDB_PREDICTION_URL = f"{_AFDB_BASE}/api/prediction"
_PUBCHEM_PROPERTY_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
_PUBCHEM_CID_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
_UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
_HUMAN_TAXONOMY_ID = 9606
_TIMEOUT_SECONDS = 30
_MAX_DOWNLOAD_BYTES = 100 * 1024 * 1024
_MAX_PUBCHEM_SDF_BYTES = 20 * 1024 * 1024
_USER_AGENT = "PyMOL-structure2text/1.0"

# UniProt accessions may include an isoform suffix (for example P04637-2).
_UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$",
    re.IGNORECASE,
)


class Af2Error(RuntimeError):
    """An error suitable for display in the PyMOL command line."""


class SmilesError(RuntimeError):
    """An error suitable for display by the SMILES commands."""


class PubChemError(RuntimeError):
    """An error suitable for display by the pubchem command."""


class PeptideError(RuntimeError):
    """An error suitable for display by the peptide command."""


class TargetError(RuntimeError):
    """An error suitable for display by the target command."""


def _request(url: str):
    request = Request(url, headers={"Accept": "application/json", "User-Agent": _USER_AGENT})
    try:
        return urlopen(request, timeout=_TIMEOUT_SECONDS)
    except HTTPError as exc:
        if exc.code == 404:
            raise Af2Error("No matching AlphaFold DB entry was found.") from exc
        if exc.code == 400:
            raise Af2Error("The search query or UniProt accession is invalid.") from exc
        raise Af2Error(f"AlphaFold DB returned HTTP {exc.code}.") from exc
    except (URLError, TimeoutError) as exc:
        raise Af2Error(f"Could not connect to AlphaFold DB: {exc}") from exc


def _request_json(url: str) -> Any:
    with _request(url) as response:
        try:
            return json.load(response)
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise Af2Error("AlphaFold DB returned an invalid JSON response.") from exc


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _has_human_taxon(record: dict[str, Any]) -> bool:
    tax_ids = _as_list(record.get("taxId")) + _as_list(record.get("taxId_display"))
    if any(str(tax_id) == str(_HUMAN_TAXONOMY_ID) for tax_id in tax_ids):
        return True
    organisms = _as_list(record.get("organismScientificName"))
    return any(str(name).casefold() == "homo sapiens" for name in organisms)


def _first(record: dict[str, Any], *keys: str) -> str:
    for key in keys:
        values = _as_list(record.get(key))
        if values:
            return str(values[0])
    return ""


def _solr_escape(value: str) -> str:
    """Escape user text before placing it in the AFDB search expression."""

    # The website search endpoint uses a Solr-style query.  Escaping every
    # non-word character prevents a PyMOL argument from changing the query.
    return re.sub(r"([^A-Za-z0-9_])", r"\\\1", value)


def _search_human_entry(query: str) -> dict[str, Any]:
    if any(character.isspace() for character in query):
        # Solr accepts a quoted phrase for multi-word protein names. Escaped
        # spaces in a wildcard expression are rejected by this endpoint.
        escaped_phrase = query.replace("\\", "\\\\").replace('"', '\\"')
        search_expression = f'(text:"{escaped_phrase}")'
    else:
        escaped_query = _solr_escape(query)
        search_expression = f"(text:*{escaped_query} OR text:{escaped_query}*)"
    params = urlencode(
        {
            "q": search_expression,
            "type": "main",
            "start": 0,
            "rows": 100,
        }
    )
    payload = _request_json(f"{_AFDB_SEARCH_URL}?{params}")
    if not isinstance(payload, dict):
        raise Af2Error("AlphaFold DB returned an unexpected search response.")

    docs = payload.get("docs", [])
    if not isinstance(docs, list):
        raise Af2Error("AlphaFold DB returned an unexpected search response.")

    candidates = [
        doc
        for doc in docs
        if isinstance(doc, dict)
        and _has_human_taxon(doc)
        and doc.get("isUniProt", True)
        and not doc.get("isComplex", False)
    ]
    if not candidates:
        raise Af2Error(f"No human protein matched '{query}'.")

    query_folded = query.casefold()

    def rank(item: tuple[int, dict[str, Any]]) -> tuple[int, ...]:
        index, doc = item
        accessions = [str(value).casefold() for value in _as_list(doc.get("uniprotAccession"))]
        genes = [str(value).casefold() for value in _as_list(doc.get("gene"))]
        descriptions = [
            str(value).casefold() for value in _as_list(doc.get("uniprotDescription"))
        ]
        accession = _first(doc, "uniprotAccession", "uniprotAccession_display")
        is_isoform = bool(doc.get("isIsoform", "-" in accession))
        exact_gene = query_folded in genes
        exact_description = query_folded in descriptions
        return (
            int(query_folded in accessions),
            int(exact_gene or exact_description),
            int(exact_gene),
            int(bool(doc.get("isUniProtReviewed", doc.get("isReviewed", False)))),
            int(
                bool(
                    doc.get(
                        "isUniProtReferenceProteome",
                        doc.get("isReferenceProteome", False),
                    )
                )
            ),
            int(not is_isoform),
            int(exact_description),
            -index,  # Preserve API relevance order when all other fields tie.
        )

    return max(enumerate(candidates), key=rank)[1]


def _prediction_records(accession: str) -> list[dict[str, Any]]:
    url = f"{_AFDB_PREDICTION_URL}/{quote(accession, safe='-')}"
    payload = _request_json(url)
    if not isinstance(payload, list):
        raise Af2Error("AlphaFold DB returned an unexpected prediction response.")
    records = [record for record in payload if isinstance(record, dict)]
    if not records:
        raise Af2Error(f"No predicted structure was found for {accession}.")
    return records


def _select_human_prediction(
    records: Iterable[dict[str, Any]], accession: str
) -> dict[str, Any]:
    accession_folded = accession.casefold()
    human_records = [record for record in records if _has_human_taxon(record)]
    if not human_records:
        raise Af2Error(f"{accession} is not a human protein.")

    exact = [
        record
        for record in human_records
        if _first(record, "uniprotAccession").casefold() == accession_folded
    ]
    if exact:
        return exact[0]

    # An accession request can return canonical and isoform records together.
    canonical = [
        record
        for record in human_records
        if "-" not in _first(record, "uniprotAccession")
    ]
    return canonical[0] if canonical else human_records[0]


def _resolve_prediction(query: str) -> dict[str, Any]:
    normalized = query.strip()
    if not normalized:
        raise Af2Error("Specify a search query. Example: af2 TP53")

    if _UNIPROT_ACCESSION_RE.fullmatch(normalized):
        accession = normalized.upper()
    else:
        search_result = _search_human_entry(normalized)
        accession = _first(
            search_result,
            "uniprotAccession",
            "uniprotAccession_display",
        )
        if not accession:
            raise Af2Error("Could not obtain a UniProt accession from the search result.")

    return _select_human_prediction(_prediction_records(accession), accession)


def _validated_pdb_url(record: dict[str, Any]) -> str:
    pdb_url = record.get("pdbUrl")
    if not isinstance(pdb_url, str) or not pdb_url:
        raise Af2Error("This entry does not provide a downloadable PDB file.")

    parsed = urlparse(pdb_url)
    if (
        parsed.scheme != "https"
        or parsed.hostname != "alphafold.ebi.ac.uk"
        or not parsed.path.startswith("/files/")
        or not parsed.path.lower().endswith(".pdb")
    ):
        raise Af2Error("AlphaFold DB returned an unsafe PDB URL.")
    return pdb_url


def _download_pdb(pdb_url: str) -> Path:
    filename = Path(urlparse(pdb_url).path).name
    if not filename or filename in {".", ".."}:
        raise Af2Error("Could not determine the PDB filename.")

    destination = Path.cwd() / filename
    if destination.exists():
        if not destination.is_file():
            raise Af2Error(f"The destination is not a regular file: {destination}")
        print(f"af2: Reusing existing file: {destination}")
        return destination

    with _request(pdb_url) as response:
        final_url = response.geturl()
        if urlparse(final_url).hostname != "alphafold.ebi.ac.uk":
            raise Af2Error("The PDB download was redirected to an unexpected host.")
        data = response.read(_MAX_DOWNLOAD_BYTES + 1)

    if len(data) > _MAX_DOWNLOAD_BYTES:
        raise Af2Error("The PDB file exceeds the 100 MiB size limit.")
    if not data:
        raise Af2Error("AlphaFold DB returned an empty PDB file.")

    try:
        with destination.open("xb") as handle:
            handle.write(data)
    except FileExistsError:
        # Another command may have completed the same download concurrently.
        print(f"af2: Reusing existing file: {destination}")
    except OSError as exc:
        raise Af2Error(f"Could not save the PDB file: {exc}") from exc

    return destination


def _numbered_object_name(prefix: str, query: str) -> str:
    target_name = re.sub(r"[^A-Za-z0-9_]+", "_", query.strip()).strip("_")
    target_name = target_name or "target"
    existing_names = set(cmd.get_names("objects"))

    for index in range(1, 100):
        candidate = f"{prefix}_{target_name}_{index:02d}"
        if candidate not in existing_names:
            return candidate

    raise Af2Error(
        f"No available object name remains in the {prefix}_{target_name}_xx series."
    )


def _object_name(query: str) -> str:
    return _numbered_object_name("AF", query)


def af2(query: str, object_name: str = "") -> Optional[str]:
    """Download and load a human AlphaFold model.

    ``query`` may be a gene name, protein name, or UniProt accession.  An
    optional second PyMOL argument sets the loaded object's name::

        af2 TP53, my_p53
    """

    try:
        record = _resolve_prediction(str(query))
        accession = _first(record, "uniprotAccession")
        gene = _first(record, "gene") or "-"
        description = _first(record, "uniprotDescription") or "-"
        pdb_path = _download_pdb(_validated_pdb_url(record))
        pymol_name = str(object_name).strip() or _object_name(str(query))
        cmd.load(str(pdb_path), pymol_name)
        print(
            f"af2: Loaded {description} (gene={gene}, UniProt={accession}) "
            f"as '{pymol_name}'."
        )
        return str(pdb_path)
    except Af2Error as exc:
        print(f"af2: Error: {exc}")
        return None


cmd.extend("af2", af2)


def _selection_to_smiles(selection: str) -> str:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise SmilesError("RDKit is not installed.") from exc

    hetatm_selection = f"({selection}) and hetatm"
    try:
        atom_count = cmd.count_atoms(hetatm_selection)
    except Exception as exc:
        raise SmilesError(f"Could not interpret selection '{selection}'.") from exc

    if atom_count == 0:
        raise SmilesError(f"Selection '{selection}' contains no HETATM records.")

    temporary_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temporary_file:
            temporary_path = Path(temporary_file.name)

        cmd.save(
            str(temporary_path),
            hetatm_selection,
            state=-1,
            format="sdf",
            quiet=1,
        )
        molecule = Chem.MolFromMolFile(
            str(temporary_path),
            sanitize=True,
            removeHs=True,
            strictParsing=False,
        )
        if molecule is None:
            raise SmilesError("Could not construct a molecule from the HETATM records.")

        smiles_text = Chem.MolToSmiles(molecule, canonical=True)
        if not smiles_text:
            raise SmilesError("Could not generate a SMILES string.")
        return smiles_text
    except SmilesError:
        raise
    except Exception as exc:
        raise SmilesError(f"SMILES conversion failed: {exc}") from exc
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass


def _copy_to_clipboard(text: str) -> None:
    if sys.platform == "darwin":
        candidates = [["pbcopy"]]
    elif sys.platform == "win32":
        candidates = [["clip.exe"]]
    else:
        candidates = [
            ["wl-copy"],
            ["xclip", "-selection", "clipboard"],
            ["xsel", "--clipboard", "--input"],
        ]

    for command in candidates:
        if shutil.which(command[0]) is None:
            continue
        try:
            subprocess.run(
                command,
                input=text,
                text=True,
                check=True,
                capture_output=True,
            )
            return
        except (OSError, subprocess.CalledProcessError):
            continue

    raise SmilesError("No supported clipboard command is available.")


def printsmiles(selection: str = "sele") -> Optional[str]:
    """Print and copy the SMILES of HETATM records in a PyMOL selection."""

    target = str(selection).strip() or "sele"
    try:
        smiles_text = _selection_to_smiles(target)
        print(f"SMILES:\n{smiles_text}")
        _copy_to_clipboard(smiles_text)
        print("printsmiles: Copied to the clipboard.")
        return smiles_text
    except SmilesError as exc:
        print(f"printsmiles: Error: {exc}")
        return None


cmd.extend("printsmiles", printsmiles)


def _smiles_object_name(requested_name: str) -> str:
    base_name = requested_name.strip() or "SMILES"
    base_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", base_name).strip("._")
    if not base_name:
        raise SmilesError("The object name is empty or invalid.")

    existing_names = {name.casefold() for name in cmd.get_names("objects")}
    if base_name.casefold() not in existing_names:
        return base_name
    for index in range(2, 100):
        candidate = f"{base_name}_{index:02d}"
        if candidate.casefold() not in existing_names:
            return candidate
    raise SmilesError(f"No available object name remains for '{base_name}'.")


def smiles(smiles_string: str = "", name: Optional[str] = None) -> Optional[str]:
    """Build and load a 3D object from a SMILES string."""

    smiles_text = str(smiles_string).strip()
    if not smiles_text:
        print("smiles: Error: Specify a SMILES string.")
        return None

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("smiles: Error: RDKit is not installed.")
        return None

    try:
        molecule = Chem.MolFromSmiles(smiles_text)
        if molecule is None:
            raise SmilesError("The supplied SMILES string is invalid.")
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        molecule = Chem.AddHs(molecule)

        embedded = False
        for random_seed in (42, 43):
            parameters = AllChem.ETKDGv3()
            parameters.randomSeed = random_seed
            parameters.useRandomCoords = True
            if AllChem.EmbedMolecule(molecule, parameters) == 0:
                embedded = True
                break
        if not embedded:
            raise SmilesError("RDKit could not generate a 3D conformation.")

        optimization_warning = ""
        if AllChem.UFFHasAllMoleculeParams(molecule):
            if AllChem.UFFOptimizeMolecule(molecule, maxIters=1000) != 0:
                optimization_warning = "UFF optimization did not fully converge."
        elif AllChem.MMFFHasAllMoleculeParams(molecule):
            if AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000) != 0:
                optimization_warning = "MMFF optimization did not fully converge."
        else:
            optimization_warning = "No supported force-field parameters were available."

        object_name = _smiles_object_name(str(name or ""))
        pdb_block = Chem.MolToPDBBlock(molecule)
        if not pdb_block.strip():
            raise SmilesError("RDKit could not generate a PDB block.")

        cmd.read_pdbstr(pdb_block, object_name)
        cmd.show("sticks", object_name)
        cmd.hide("everything", f"({object_name} and hydro and (elem C extend 1))")
        cmd.orient(object_name)
        print(f"SMILES:\n{canonical_smiles}")
        print(f"smiles: Loaded the 3D object as '{object_name}'.")
        if optimization_warning:
            print(f"smiles: Warning: {optimization_warning}")
        return canonical_smiles
    except SmilesError as exc:
        print(f"smiles: Error: {exc}")
        return None
    except Exception as exc:
        print(f"smiles: Error: Could not build the 3D object: {exc}")
        return None


cmd.extend("smiles", smiles)


def _pubchem_smiles(query: str):
    encoded_query = quote(query, safe="")
    url = f"{_PUBCHEM_PROPERTY_URL}/{encoded_query}/property/SMILES,ConnectivitySMILES/JSON"
    request = Request(url, headers={"Accept": "application/json", "User-Agent": _USER_AGENT})
    try:
        with urlopen(request, timeout=_TIMEOUT_SECONDS) as response:
            payload = json.load(response)
    except HTTPError as exc:
        if exc.code == 404:
            raise PubChemError(f"No PubChem compound matched '{query}'.") from exc
        raise PubChemError(f"PubChem returned HTTP {exc.code}.") from exc
    except (URLError, TimeoutError) as exc:
        raise PubChemError(f"Could not connect to PubChem: {exc}") from exc
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise PubChemError("PubChem returned an invalid JSON response.") from exc

    try:
        properties = payload["PropertyTable"]["Properties"]
        first_result = properties[0]
    except (KeyError, IndexError, TypeError) as exc:
        raise PubChemError("PubChem returned an unexpected property response.") from exc

    smiles_text = first_result.get("SMILES") or first_result.get("ConnectivitySMILES")
    if not smiles_text:
        raise PubChemError(f"PubChem returned no SMILES for '{query}'.")

    return first_result.get("CID"), str(smiles_text)


def _pubchem_3d_sdf(cid: Any) -> bytes:
    if cid is None:
        raise PubChemError("PubChem returned no CID for the selected compound.")

    url = f"{_PUBCHEM_CID_URL}/{quote(str(cid), safe='')}/SDF?record_type=3d"
    request = Request(
        url,
        headers={"Accept": "chemical/x-mdl-sdfile", "User-Agent": _USER_AGENT},
    )
    try:
        with urlopen(request, timeout=_TIMEOUT_SECONDS) as response:
            if urlparse(response.geturl()).hostname != "pubchem.ncbi.nlm.nih.gov":
                raise PubChemError("The PubChem SDF download was redirected unexpectedly.")
            sdf_data = response.read(_MAX_PUBCHEM_SDF_BYTES + 1)
    except HTTPError as exc:
        if exc.code in {400, 404}:
            raise PubChemError(f"No PubChem 3D conformer is available for CID {cid}.") from exc
        raise PubChemError(f"PubChem returned HTTP {exc.code} for the 3D SDF.") from exc
    except (URLError, TimeoutError) as exc:
        raise PubChemError(f"Could not download the PubChem 3D SDF: {exc}") from exc

    if len(sdf_data) > _MAX_PUBCHEM_SDF_BYTES:
        raise PubChemError("The PubChem 3D SDF exceeds the 20 MiB size limit.")
    if not sdf_data or b"$$$$" not in sdf_data:
        raise PubChemError("PubChem returned an invalid 3D SDF.")
    return sdf_data


def _load_pubchem_sdf(sdf_data: bytes, object_name: str) -> None:
    temporary_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temporary_file:
            temporary_file.write(sdf_data)
            temporary_path = Path(temporary_file.name)
        cmd.load(str(temporary_path), object_name)
    except OSError as exc:
        raise PubChemError(f"Could not prepare the PubChem 3D SDF: {exc}") from exc
    except Exception as exc:
        raise PubChemError(f"Could not load the PubChem 3D SDF into PyMOL: {exc}") from exc
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass


def pubchem(query: str = "") -> Optional[str]:
    """Search PubChem, copy the SMILES, and load its 3D SDF conformer."""

    query_text = str(query).strip()
    if not query_text:
        print("pubchem: Error: Specify a PubChem search query.")
        return None

    try:
        cid, smiles_text = _pubchem_smiles(query_text)
        sdf_data = _pubchem_3d_sdf(cid)
        object_name = _numbered_object_name("PC", query_text)
        _load_pubchem_sdf(sdf_data, object_name)
        cmd.hide("everything", f"({object_name} and hydro and (elem C extend 1))")

        print(f"SMILES:\n{smiles_text}")
        try:
            _copy_to_clipboard(smiles_text)
            print("pubchem: Copied the SMILES to the clipboard.")
        except SmilesError as exc:
            print(
                "pubchem: Warning: Could not copy the SMILES to the "
                f"clipboard: {exc}"
            )

        cid_text = str(cid) if cid is not None else "unknown"
        print(
            f"pubchem: Loaded the PubChem 3D SDF for CID {cid_text} "
            f"as '{object_name}'."
        )
        return smiles_text
    except PubChemError as exc:
        print(f"pubchem: Error: {exc}")
        return None
    except Exception as exc:
        print(f"pubchem: Error: Could not load the structure into PyMOL: {exc}")
        return None


cmd.extend("pubchem", pubchem)


_PEPTIDE_SIDE_CHAINS = {
    "G": "([H])([H])",
    "A": "([H])(C)",
    "V": "([H])(C(C)C)",
    "L": "([H])(CC(C)C)",
    "I": "([H])(C(C)CC)",
    "M": "([H])(CCSC)",
    "F": "([H])(CC1=CC=CC=C1)",
    "W": "([H])(CC1=CNC2=C1C=CC=C2)",
    "S": "([H])(CO)",
    "C": "([H])(CS)",
    "Y": "([H])(CC1=CC=C(C=C1)O)",
    "N": "([H])(C(C(=O)N))",
    "Q": "([H])(C(CC(=O)N))",
    "D": "([H])(C(C(=O)O))",
    "E": "([H])(C(CC(=O)O))",
    "H": "([H])(C(C1=CNC=N1))",
    "K": "([H])(CCCCN)",
    "R": "([H])(CCCNC(=N)N)",
    "T": "([H])([C@]([H])(O)C)",
    "a": "(C)([H])",
    "v": "(C(C)C)([H])",
    "l": "(CC(C)C)([H])",
    "i": "(C(C)CC)([H])",
    "m": "(CCSC)([H])",
    "f": "(CC1=CC=CC=C1)([H])",
    "w": "(CC1=CNC2=C1C=CC=C2)([H])",
    "s": "(CO)([H])",
    "c": "(CS)([H])",
    "y": "(CC1=CC=C(C=C1)O)([H])",
    "n": "(C(C(=O)N))([H])",
    "q": "(C(CC(=O)N))([H])",
    "d": "(C(C(=O)O))([H])",
    "e": "(C(CC(=O)O))([H])",
    "h": "(C(C1=CNC=N1))([H])",
    "k": "(CCCCN)([H])",
    "r": "(CCCNC(=N)N)([H])",
    "t": "([C@]([H])(O)C)([H])",
    "g": "(C)(C)",
    "J": "([H])(CS8)",
    "j": "(CS8)([H])",
    "B": "([H])(CS7)",
    "b": "(CS7)([H])",
}
_PEPTIDE_RESIDUES = {
    code: f"N[C@@]{side_chain}C(=O)"
    for code, side_chain in _PEPTIDE_SIDE_CHAINS.items()
}
_PEPTIDE_RESIDUES.update(
    {
        "P": "N1[C@@]([H])(CCC1)C(=O)",
        "p": "N1[C@@](CCC1)([H])C(=O)",
    }
)
_PEPTIDE_MODES = {
    "0": "linear",
    "1": "head-to-tail",
    "2": "simple disulfide",
    "3": "cysteine cyclization",
}


def _peptide_smiles(sequence: str, mode: str) -> str:
    invalid_codes = sorted(set(sequence) - set(_PEPTIDE_RESIDUES))
    if invalid_codes:
        invalid_text = ", ".join(repr(code) for code in invalid_codes)
        raise PeptideError(f"Unsupported residue code(s): {invalid_text}.")

    for left, right in (("B", "b"), ("J", "j")):
        ring_marker_count = sequence.count(left) + sequence.count(right)
        if ring_marker_count % 2:
            raise PeptideError(
                f"Residue codes {left}/{right} must occur in pairs."
            )

    residue_smiles = "".join(_PEPTIDE_RESIDUES[code] for code in sequence)
    if mode == "1":
        return "N9" + residue_smiles[1:-4] + "9(=O)"
    if mode == "2":
        return (
            "N[C@@]([H])(CS9)C(=O)"
            + residue_smiles
            + "N[C@@]([H])(CS9)C(=O)O"
        )
    if mode == "3":
        return (
            "C9C(=O)"
            + residue_smiles
            + "N[C@@]([H])(CS9)C(=O)N"
        )
    return residue_smiles + "N"


def _peptide_object_name(sequence: str, requested_name: str) -> str:
    base_name = requested_name.strip() or sequence
    base_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", base_name).strip("._")
    if not base_name:
        raise PeptideError("The object name is empty or invalid.")

    existing_names = {name.casefold() for name in cmd.get_names("objects")}
    if base_name.casefold() not in existing_names:
        return base_name
    for index in range(2, 100):
        candidate = f"{base_name}_{index:02d}"
        if candidate.casefold() not in existing_names:
            return candidate
    raise PeptideError(f"No available object name remains for '{base_name}'.")


def peptide(
    sequence: str = "",
    mode: Optional[str] = None,
    name: Optional[str] = None,
) -> Optional[str]:
    """Build and load a 3D peptide from one-letter residue codes."""

    sequence_text = str(sequence).strip()
    if not sequence_text:
        print("peptide: Error: Specify a peptide sequence.")
        return None

    mode_text = "0" if mode is None or not str(mode).strip() else str(mode).strip()
    if mode_text not in _PEPTIDE_MODES:
        print("peptide: Error: Mode must be 0, 1, 2, or 3.")
        return None

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("peptide: Error: RDKit is not installed.")
        return None

    try:
        smiles_text = _peptide_smiles(sequence_text, mode_text)
        molecule = Chem.MolFromSmiles(smiles_text)
        if molecule is None:
            raise PeptideError("The generated peptide SMILES is invalid.")
        molecule = Chem.AddHs(molecule)

        embedded = False
        for random_seed in (42, 43):
            parameters = AllChem.ETKDGv3()
            parameters.randomSeed = random_seed
            parameters.useRandomCoords = True
            if AllChem.EmbedMolecule(molecule, parameters) == 0:
                embedded = True
                break
        if not embedded:
            raise PeptideError("RDKit could not generate a 3D peptide conformation.")

        optimization_status = AllChem.UFFOptimizeMolecule(molecule, maxIters=1000)
        object_name = _peptide_object_name(sequence_text, str(name or ""))
        pdb_block = Chem.MolToPDBBlock(molecule)
        if not pdb_block.strip():
            raise PeptideError("RDKit could not generate a PDB block.")

        cmd.read_pdbstr(pdb_block, object_name)
        cmd.show("sticks", object_name)
        cmd.orient(object_name)
        print(f"SMILES:\n{smiles_text}")
        print(
            f"peptide: Loaded {_PEPTIDE_MODES[mode_text]} peptide "
            f"'{sequence_text}' as '{object_name}'."
        )
        if optimization_status != 0:
            print("peptide: Warning: UFF optimization did not fully converge.")
        return smiles_text
    except PeptideError as exc:
        print(f"peptide: Error: {exc}")
        return None
    except Exception as exc:
        print(f"peptide: Error: Could not build the peptide: {exc}")
        return None


cmd.extend("peptide", peptide)


def mkdir(directory: str = "") -> Optional[str]:
    """Create a directory, including any missing parent directories."""

    directory_text = str(directory).strip()
    if not directory_text:
        print("mkdir: Error: Specify a directory path.")
        return None

    path = Path(directory_text)
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"mkdir: Error: Path already exists: {path}")
        return None
    except OSError as exc:
        print(f"mkdir: Error: Could not create '{path}': {exc}")
        return None

    resolved_path = path.resolve()
    print(f"mkdir: Created directory: {resolved_path}")
    return str(resolved_path)


cmd.extend("mkdir", mkdir)


def rm(filename: str = "") -> Optional[Any]:
    """Delete one file or files matching a non-recursive glob pattern."""

    filename_text = str(filename).strip()
    if not filename_text:
        print("rm: Error: Specify a file path.")
        return None
    if "**" in filename_text:
        print("rm: Error: Recursive glob patterns are not supported.")
        return None

    if glob.has_magic(filename_text):
        matches = sorted(glob.glob(filename_text))
        if not matches:
            print(f"rm: Error: No files matched pattern: {filename_text}")
            return None
    else:
        matches = [filename_text]

    deleted_paths = []
    for match in matches:
        path = Path(match)
        if not path.exists() and not path.is_symlink():
            print(f"rm: Error: File does not exist: {path}")
            continue
        if path.is_dir() and not path.is_symlink():
            print(f"rm: Error: Refusing to delete a directory: {path}")
            continue

        resolved_path = path.absolute()
        try:
            path.unlink()
        except OSError as exc:
            print(f"rm: Error: Could not delete '{path}': {exc}")
            continue

        print(f"rm: Deleted file: {resolved_path}")
        deleted_paths.append(str(resolved_path))

    if not deleted_paths:
        return None
    if len(deleted_paths) == 1:
        return deleted_paths[0]
    return deleted_paths


cmd.extend("rm", rm)


def _uniprot_query_json(query: str) -> dict[str, Any]:
    params = urlencode({"query": query, "format": "json", "size": 25})
    request = Request(
        f"{_UNIPROT_SEARCH_URL}?{params}",
        headers={"Accept": "application/json", "User-Agent": _USER_AGENT},
    )
    try:
        with urlopen(request, timeout=_TIMEOUT_SECONDS) as response:
            payload = json.load(response)
    except HTTPError as exc:
        raise TargetError(f"UniProt returned HTTP {exc.code}.") from exc
    except (URLError, TimeoutError) as exc:
        raise TargetError(f"Could not connect to UniProt: {exc}") from exc
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise TargetError("UniProt returned an invalid JSON response.") from exc

    if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
        raise TargetError("UniProt returned an unexpected search response.")
    return payload


def _uniprot_gene_names(record: dict[str, Any]) -> list[str]:
    names = []
    for gene in record.get("genes", []):
        if not isinstance(gene, dict):
            continue
        gene_name = gene.get("geneName", {})
        if isinstance(gene_name, dict) and gene_name.get("value"):
            names.append(str(gene_name["value"]))
        for synonym in gene.get("synonyms", []):
            if isinstance(synonym, dict) and synonym.get("value"):
                names.append(str(synonym["value"]))
    return names


def _uniprot_target_record(query: str) -> dict[str, Any]:
    escaped = query.replace("\\", "\\\\").replace('"', '\\"')
    if _UNIPROT_ACCESSION_RE.fullmatch(query):
        queries = [f"(accession:{query.upper()}) AND (organism_id:9606)"]
    else:
        queries = [
            f'(gene_exact:"{escaped}") AND (organism_id:9606) AND (reviewed:true)',
            f'(gene_exact:"{escaped}") AND (organism_id:9606)',
            f'("{escaped}") AND (organism_id:9606) AND (reviewed:true)',
            f'("{escaped}") AND (organism_id:9606)',
        ]

    results = []
    for search_query in queries:
        results = _uniprot_query_json(search_query)["results"]
        if results:
            break
    if not results:
        raise TargetError(f"No human UniProt entry matched '{query}'.")

    query_folded = query.casefold()

    def rank(item: tuple[int, dict[str, Any]]) -> tuple[int, int, float, int]:
        index, record = item
        genes = [name.casefold() for name in _uniprot_gene_names(record)]
        accession = str(record.get("primaryAccession", "")).casefold()
        reviewed = record.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
        annotation_score = float(record.get("annotationScore", 0.0) or 0.0)
        return (
            int(query_folded == accession or query_folded in genes),
            int(reviewed),
            annotation_score,
            -index,
        )

    selected = max(enumerate(results), key=rank)[1]
    organism = selected.get("organism", {})
    if not isinstance(organism, dict) or organism.get("taxonId") != 9606:
        raise TargetError("The selected UniProt entry is not a human protein.")
    return selected


def _pdb_references(record: dict[str, Any]) -> list[dict[str, str]]:
    references = []
    seen = set()
    for cross_reference in record.get("uniProtKBCrossReferences", []):
        if not isinstance(cross_reference, dict) or cross_reference.get("database") != "PDB":
            continue
        pdb_id = str(cross_reference.get("id", "")).upper()
        if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id) or pdb_id in seen:
            continue
        properties = {
            str(item.get("key")): str(item.get("value"))
            for item in cross_reference.get("properties", [])
            if isinstance(item, dict) and item.get("key") and item.get("value")
        }
        references.append(
            {
                "pdb_id": pdb_id,
                "method": properties.get("Method", ""),
                "resolution": properties.get("Resolution", ""),
                "chains": properties.get("Chains", ""),
            }
        )
        seen.add(pdb_id)
    return references


def _write_structure_list(
    path: Path,
    accession: str,
    references: list[dict[str, str]],
) -> None:
    lines = ["uniprot_accession\tpdb_id\tmethod\tresolution\tchains"]
    for reference in references:
        values = [
            accession,
            reference["pdb_id"],
            reference["method"],
            reference["resolution"],
            reference["chains"],
        ]
        lines.append("\t".join(value.replace("\t", " ").replace("\n", " ") for value in values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _largest_hetatm_component(object_name: str) -> Optional[dict[str, Any]]:
    model = cmd.get_model(f"{object_name} and hetatm and not solvent", state=1)
    atoms = list(getattr(model, "atom", []))
    if not atoms:
        return None

    parents = list(range(len(atoms)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    for bond in getattr(model, "bond", []):
        indices = getattr(bond, "index", ())
        if len(indices) == 2:
            left, right = int(indices[0]), int(indices[1])
            if 0 <= left < len(atoms) and 0 <= right < len(atoms):
                union(left, right)

    components: dict[int, list[Any]] = {}
    for index, atom in enumerate(atoms):
        components.setdefault(find(index), []).append(atom)

    def heavy_atom_count(component: list[Any]) -> int:
        count = 0
        for atom in component:
            element = getattr(atom, "symbol", "") or getattr(atom, "elem", "")
            if str(element).upper() != "H":
                count += 1
        return count

    largest = max(components.values(), key=lambda component: (heavy_atom_count(component), len(component)))
    heavy_count = heavy_atom_count(largest)
    if heavy_count == 0:
        return None

    chains = Counter(str(getattr(atom, "chain", "")) for atom in largest)
    residue_names = sorted({str(getattr(atom, "resn", "")) for atom in largest})
    residue_ids = sorted({str(getattr(atom, "resi", "")) for atom in largest})
    return {
        "atom_indices": [int(atom.index) for atom in largest],
        "chain": chains.most_common(1)[0][0],
        "residue_names": ",".join(name for name in residue_names if name),
        "residue_ids": ",".join(value for value in residue_ids if value),
        "heavy_atom_count": heavy_count,
    }


def _select_ligand_chain(
    object_name: str,
    ligand_selection: str,
    ligand_chain: str,
) -> tuple[str, str, int]:
    protein_model = cmd.get_model(
        f"{object_name} and polymer.protein and name CA",
        state=1,
    )
    protein_atoms = list(getattr(protein_model, "atom", []))
    if not protein_atoms:
        raise TargetError(f"{object_name} contains no protein CA atoms.")

    chain_counts = Counter(str(getattr(atom, "chain", "")) for atom in protein_atoms)
    selected_chain = ligand_chain if chain_counts.get(ligand_chain, 0) >= 3 else ""

    if not selected_chain:
        contact_model = cmd.get_model(
            f"({object_name} and polymer.protein and name CA) within 8 of {ligand_selection}",
            state=1,
        )
        contact_counts = Counter(
            str(getattr(atom, "chain", "")) for atom in getattr(contact_model, "atom", [])
        )
        if contact_counts:
            selected_chain = contact_counts.most_common(1)[0][0]
        else:
            selected_chain = chain_counts.most_common(1)[0][0]

    selected_indices = [
        int(atom.index)
        for atom in protein_atoms
        if str(getattr(atom, "chain", "")) == selected_chain
    ]
    if len(selected_indices) < 3:
        raise TargetError(f"{object_name} has too few CA atoms for alignment.")

    selection_name = f"align_{object_name}"
    index_expression = "+".join(str(index) for index in selected_indices)
    cmd.select(selection_name, f"{object_name} and index {index_expression}")
    return selection_name, selected_chain, len(selected_indices)


def _target_safe_name(query: str) -> str:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", query.strip()).strip("._")
    if not safe_name or safe_name in {".", ".."}:
        raise TargetError("The target name cannot be converted to a safe directory name.")
    return safe_name


def target(query: str = "") -> Optional[str]:
    """Fetch and align ligand-bound PDB structures for a human UniProt target."""

    query_text = str(query).strip()
    if not query_text:
        print("target: Error: Specify a target name.")
        return None

    try:
        safe_name = _target_safe_name(query_text)
    except TargetError as exc:
        print(f"target: Error: {exc}")
        return None

    target_directory = Path.cwd() / safe_name
    if target_directory.exists():
        print(f"target: Error: Directory already exists: {target_directory}")
        return None

    try:
        target_directory.mkdir(parents=False, exist_ok=False)
        os.chdir(target_directory)
    except OSError as exc:
        print(f"target: Error: Could not create or enter '{target_directory}': {exc}")
        return None

    log_path = target_directory / f"{safe_name}.log"
    session_path = target_directory / f"{safe_name}.pse"
    structure_list_path = target_directory / "structures.tsv"
    completed = False
    temporary_selections = []

    with log_path.open("w", encoding="utf-8") as log_handle:

        def log(message: str, level: str = "INFO") -> None:
            timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
            log_handle.write(f"{timestamp}\t{level}\t{message}\n")
            log_handle.flush()
            print(f"target: {message}")

        log(f"Created and entered directory: {target_directory}")
        try:
            record = _uniprot_target_record(query_text)
            accession = str(record.get("primaryAccession", ""))
            uniprot_id = str(record.get("uniProtkbId", ""))
            genes = _uniprot_gene_names(record)
            log(
                f"Selected UniProt entry {accession} ({uniprot_id}); "
                f"genes={','.join(genes) or 'unknown'}."
            )

            references = _pdb_references(record)
            _write_structure_list(structure_list_path, accession, references)
            log(f"Saved {len(references)} PDB cross-references to {structure_list_path.name}.")
            if not references:
                raise TargetError(f"UniProt entry {accession} has no PDB cross-references.")

            combined_site = f"target_{safe_name}_site"
            cmd.select(combined_site, "none")
            temporary_selections.append(combined_site)
            processed = []

            for reference in references:
                pdb_id = reference["pdb_id"]
                object_name = pdb_id
                try:
                    existing_objects = {
                        name.casefold() for name in cmd.get_names("objects")
                    }
                    if object_name.casefold() in existing_objects:
                        log(
                            f"Skipped {pdb_id}: object name already exists.",
                            "WARNING",
                        )
                        continue
                    cmd.fetch(
                        pdb_id.lower(),
                        object_name,
                        type="cif",
                        async_=0,
                        path=str(target_directory),
                        quiet=1,
                    )
                    if cmd.count_atoms(object_name) == 0:
                        raise TargetError("the fetched object contains no atoms")
                    log(f"Fetched {pdb_id} as {object_name}.")

                    ligand = _largest_hetatm_component(object_name)
                    if ligand is None:
                        log(f"Skipped {pdb_id}: no non-solvent HETATM molecule.", "WARNING")
                        continue

                    ligand_selection = f"lig_{object_name}"
                    ligand_indices = "+".join(str(index) for index in ligand["atom_indices"])
                    cmd.select(
                        ligand_selection,
                        f"{object_name} and index {ligand_indices}",
                    )
                    chain_selection, chain_id, ca_count = _select_ligand_chain(
                        object_name,
                        ligand_selection,
                        str(ligand["chain"]),
                    )
                    temporary_selections.append(chain_selection)
                    site_selection = f"site_{object_name}"
                    cmd.select(
                        site_selection,
                        f"({ligand_selection}) or byres "
                        f"(({object_name} and polymer.protein) within 8 of {ligand_selection})",
                    )
                    temporary_selections.append(site_selection)
                    cmd.select(combined_site, f"{combined_site} or {site_selection}")
                    cmd.show("sticks", ligand_selection)
                    cmd.show("cartoon", f"{object_name} and polymer.protein")
                    processed.append(
                        {
                            "pdb_id": pdb_id,
                            "object_name": object_name,
                            "ligand_selection": ligand_selection,
                            "chain_selection": chain_selection,
                            "chain_id": chain_id,
                        }
                    )
                    log(
                        f"Selected largest HETATM molecule in {pdb_id}: "
                        f"resn={ligand['residue_names'] or 'unknown'}, "
                        f"resi={ligand['residue_ids'] or 'unknown'}, "
                        f"heavy_atoms={ligand['heavy_atom_count']}, "
                        f"alignment_chain={chain_id or '<blank>'}, CA_atoms={ca_count}."
                    )
                except Exception as exc:
                    log(f"Skipped {pdb_id}: {exc}", "WARNING")

            if not processed:
                raise TargetError("No fetched structure contained a usable HETATM binding site.")

            reference = processed[0]
            log(
                f"Using {reference['pdb_id']} chain "
                f"{reference['chain_id'] or '<blank>'} as the alignment reference."
            )
            for item in processed[1:]:
                try:
                    result = cmd.align(
                        item["chain_selection"],
                        reference["chain_selection"],
                        cycles=5,
                        transform=1,
                        quiet=1,
                    )
                    rmsd = result[0] if isinstance(result, (tuple, list)) and result else result
                    log(
                        f"Aligned {item['pdb_id']} chain {item['chain_id'] or '<blank>'} "
                        f"to {reference['pdb_id']}; RMSD={rmsd}."
                    )
                except Exception as exc:
                    log(f"Alignment failed for {item['pdb_id']}: {exc}", "WARNING")

            cmd.zoom(combined_site, buffer=2.0)
            log(f"Zoomed to all selected ligands and protein residues within 8 A.")
            completed = True
        except TargetError as exc:
            log(f"Error: {exc}", "ERROR")
        except Exception as exc:
            log(f"Unexpected error: {exc}", "ERROR")
        finally:
            deleted_count = 0
            for selection_name in dict.fromkeys(temporary_selections):
                try:
                    cmd.delete(selection_name)
                    deleted_count += 1
                except Exception as exc:
                    log(
                        f"Could not delete selection {selection_name}: {exc}",
                        "WARNING",
                    )
            log(f"Deleted {deleted_count} alignment/site selections.")

            disabled_count = 0
            try:
                remaining_selections = list(cmd.get_names("selections"))
            except Exception as exc:
                log(f"Could not list remaining selections: {exc}", "WARNING")
                remaining_selections = []

            selection_group = "select"
            grouped_count = 0
            try:
                cmd.group(selection_group)
                for selection_name in remaining_selections:
                    cmd.group(
                        selection_group,
                        selection_name,
                        action="add",
                        quiet=1,
                    )
                    grouped_count += 1
                log(
                    f"Grouped {grouped_count} remaining selections under "
                    f"{selection_group}."
                )
            except Exception as exc:
                log(
                    f"Could not group remaining selections under "
                    f"{selection_group}: {exc}",
                    "WARNING",
                )

            for selection_name in remaining_selections:
                try:
                    cmd.disable(selection_name)
                    disabled_count += 1
                except Exception as exc:
                    log(
                        f"Could not disable selection {selection_name}: {exc}",
                        "WARNING",
                    )
            log(f"Disabled {disabled_count} remaining selections.")

            try:
                cmd.save(str(session_path))
                log(f"Saved PyMOL session: {session_path.name}")
            except Exception as exc:
                log(f"Could not save PyMOL session: {exc}", "ERROR")
                completed = False

            deleted_cif_count = 0
            cif_cleanup_failed = False
            cif_paths = []
            for pattern in ("*.cif", "*.cif.gz"):
                cif_paths.extend(target_directory.glob(pattern))
            for cif_path in dict.fromkeys(cif_paths):
                try:
                    if cif_path.is_file() or cif_path.is_symlink():
                        cif_path.unlink()
                        deleted_cif_count += 1
                except OSError as exc:
                    cif_cleanup_failed = True
                    log(f"Could not delete fetched CIF file {cif_path.name}: {exc}", "ERROR")
            log(f"Deleted {deleted_cif_count} fetched CIF files.")
            if cif_cleanup_failed:
                completed = False

    if not completed:
        print(f"target: Error: Workflow did not complete. See {log_path}")
        return None
    print(f"target: Completed successfully. Session: {session_path}")
    return str(session_path)


cmd.extend("target", target)


def _reinvent_script_path() -> Path:
    return Path(_reinvent_script_path.__code__.co_filename).resolve()


def _reinvent_python_candidates() -> list[Path]:
    candidates = []

    def add(path: Any) -> None:
        if not path:
            return
        candidate = Path(path).expanduser()
        if candidate.is_file() and candidate not in candidates:
            candidates.append(candidate)

    add(os.environ.get("REINVENT_PYTHON"))
    add(sys.executable)
    for key in ("CONDA_PREFIX", "VIRTUAL_ENV"):
        prefix = os.environ.get(key)
        if prefix:
            executable = "Scripts/python.exe" if os.name == "nt" else "bin/python"
            add(Path(prefix) / executable)

    reinvent_executable = shutil.which("reinvent")
    if reinvent_executable:
        try:
            first_line = (
                Path(reinvent_executable)
                .read_bytes()
                .splitlines()[0]
                .decode("utf-8")
            )
            if first_line.startswith("#!") and not first_line.endswith("/env python"):
                add(first_line[2:].strip().split()[0])
        except (OSError, UnicodeError, IndexError):
            pass
        python_name = "python.exe" if os.name == "nt" else "python"
        add(Path(reinvent_executable).parent / python_name)

    home = Path.home()
    python_path = Path("Scripts/python.exe") if os.name == "nt" else Path("bin/python")
    for root in (
        home / "miniforge3",
        home / "miniconda3",
        home / "anaconda3",
        home / "mambaforge",
        home / "micromamba",
        home / ".conda",
    ):
        environment_root = root / "envs"
        environments = sorted(environment_root.glob("*")) if environment_root.is_dir() else []
        for directory in (root, *environments):
            add(directory / python_path)
    return candidates


def _reinvent_runtime() -> dict[str, str]:
    check = (
        "import pathlib,torch,rdkit,reinvent.Reinvent as entry; "
        "print(pathlib.Path(entry.__file__).resolve())"
    )
    errors = []
    for python_path in _reinvent_python_candidates():
        try:
            with tempfile.TemporaryDirectory(prefix="reinvent_probe_") as probe_dir:
                result = subprocess.run(
                    [str(python_path), "-c", check],
                    cwd=probe_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=30,
                    check=False,
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            errors.append(f"{python_path}: {exc}")
            continue
        if result.returncode != 0:
            errors.append(f"{python_path}: {result.stderr.strip()[-180:]}")
            continue
        origin = Path(result.stdout.strip().splitlines()[-1])
        if not origin.is_file() or origin.name != "Reinvent.py":
            errors.append(f"{python_path}: invalid REINVENT4 module path: {origin}")
            continue
        executable_name = "reinvent.exe" if os.name == "nt" else "reinvent"
        executable = python_path.parent / executable_name
        if executable.is_file():
            return {
                "python": str(python_path),
                "executable": str(executable),
                "package": str(origin),
            }
        errors.append(f"{python_path}: REINVENT console command missing: {executable}")
    raise RuntimeError(
        "A working REINVENT4 Python environment was not found. Install "
        "REINVENT4 with its dependencies, or set REINVENT_PYTHON to that "
        "environment's Python executable. Checked:\n"
        + "\n".join(errors[:12])
    )


def _reinvent_prior_roots(runtime: dict[str, str]) -> list[Path]:
    home = Path.home()
    package_directory = Path(runtime["package"]).parent
    bases = [
        Path(value).expanduser()
        for key in ("REINVENT_PRIORS_DIR", "REINVENT4_HOME")
        for value in os.environ.get(key, "").split(os.pathsep)
        if value
    ]
    bases += [
        package_directory.parent,
        package_directory,
        _reinvent_script_path().parent,
        Path.cwd(),
        Path.cwd().parent,
        Path(runtime["python"]).parent.parent,
        home / "REINVENT4",
        home / "reinvent4",
        home / "Documents" / "REINVENT4",
        home / "Documents" / "modeling" / "REINVENT4",
        home / "apps" / "REINVENT4",
        home / "projects" / "REINVENT4",
    ]
    roots = []
    for base in bases:
        for root in (base, base / "priors", base / "models", base / "checkpoints"):
            try:
                resolved = root.resolve()
                if resolved.is_dir() and resolved not in roots:
                    roots.append(resolved)
            except OSError:
                pass
    return roots


def _reinvent_find_priors(runtime: dict[str, str]) -> dict[str, Path]:
    files = {}
    for root in _reinvent_prior_roots(runtime):
        search_directories = [root]
        search_directories += [path for path in root.glob("*/priors") if path.is_dir()]
        search_directories += [path for path in root.glob("*/models") if path.is_dir()]
        for directory in search_directories:
            for path in sorted(directory.glob("*.prior")):
                if path.is_file() and "mol2mol" in path.stem.lower():
                    files.setdefault(path.stem.lower(), path.resolve())
    return files


def reinvent_priors() -> dict[str, Path]:
    """List Mol2Mol prior files discovered for REINVENT4."""

    try:
        runtime = _reinvent_runtime()
        print(f"[REINVENT] Environment: {runtime['python']}")
        print(f"[REINVENT] Package: {runtime['package']}")
        print(f"[REINVENT] Loaded script: {_reinvent_script_path()}")
        print("[REINVENT] Searched directories:")
        for root in _reinvent_prior_roots(runtime):
            print(f"  {root}")
        priors = _reinvent_find_priors(runtime)
        if priors:
            print("[REINVENT] Found Mol2Mol priors:")
            for prior_name, prior_path in priors.items():
                print(f"  {prior_name}: {prior_path}")
        else:
            print("[REINVENT] No Mol2Mol .prior files found.")
            print("Set REINVENT_PRIORS_DIR to the directory containing prior files.")
        return priors
    except Exception as exc:
        print(f"reinvent_priors: Error: {exc}")
        return {}


def _reinvent_select_priors(value: str, runtime: dict[str, str]) -> list[Path]:
    all_priors = _reinvent_find_priors(runtime)
    if not all_priors:
        raise RuntimeError(
            "No Mol2Mol priors were found. Run reinvent_priors to inspect search paths."
        )
    if not value or value.lower() in ("all", "*"):
        return list(all_priors.values())
    chosen = []
    for name in (item.strip() for item in value.split("+") if item.strip()):
        key = Path(name).stem.lower()
        matched = all_priors.get(key)
        if matched is None:
            candidates = [path for prior_name, path in all_priors.items() if key in prior_name]
            if len(candidates) == 1:
                matched = candidates[0]
        if matched is None:
            raise ValueError(
                f"Prior {name!r} was not uniquely found. Available: "
                f"{', '.join(all_priors)}"
            )
        if matched not in chosen:
            chosen.append(matched)
    return chosen


def _reinvent_dispatch(
    selection: str = "sele",
    prior: str = "all",
    num_smiles: Any = 10000,
    temperature: Any = 1.0,
) -> str:
    number = int(num_smiles)
    temperature_value = float(temperature)
    if number < 1 or not 0 < temperature_value <= 10:
        raise ValueError(
            "num_smiles must be positive; temperature must be greater than 0 "
            "and no greater than 10."
        )
    runtime = _reinvent_runtime()
    chosen_priors = _reinvent_select_priors(str(prior), runtime)
    atoms = cmd.get_model(f"({selection})").atom
    if not atoms:
        raise ValueError(f"No atoms are present in selection: {selection}")
    residues = {
        (atom.model, atom.segi, atom.chain, atom.resi, atom.resn)
        for atom in atoms
    }
    if len(residues) != 1:
        raise ValueError("Select exactly one ligand residue from one object.")
    model, _segi, _chain, _resi, residue_name = next(iter(residues))
    match = re.search(r"(?<![a-z0-9])([0-9][a-z0-9]{3})(?![a-z0-9])", model.lower())
    pdb_id = match.group(1).upper() if match else _reinvent_safe_filename(model)
    base_name = "_".join(
        (
            dt.datetime.now().strftime("%Y%m%d"),
            _reinvent_safe_filename(pdb_id),
            _reinvent_safe_filename(residue_name),
        )
    )
    output_directory = Path.cwd() / "reinvent"
    output_directory.mkdir(parents=True, exist_ok=True)
    output_csv = output_directory / f"{base_name}.csv"
    job_directory = Path(
        tempfile.mkdtemp(prefix=f"{base_name}_", dir=output_directory)
    )
    molecule_path = job_directory / "selected_ligand.mol"
    cmd.save(str(molecule_path), selection, state=-1, format="mol")
    if not molecule_path.exists() or molecule_path.stat().st_size == 0:
        raise RuntimeError("PyMOL failed to save the selected molecule.")
    manifest_path = job_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "device": "auto",
                "priors": [str(path) for path in chosen_priors],
                "num_smiles": number,
                "temperature": temperature_value,
                "output_csv": str(output_csv),
                "runtime": runtime,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    log_path = job_directory / "job.log"
    with log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            [
                runtime["python"],
                str(_reinvent_script_path()),
                "--reinvent-worker",
                str(manifest_path),
            ],
            cwd=str(job_directory),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    print(
        f"[REINVENT] Started PID={process.pid}; priors={len(chosen_priors)}; "
        f"requested={number} per prior"
    )
    print(f"[REINVENT] Output CSV: {output_csv}")
    print(f"[REINVENT] Logs: {log_path}")
    return str(output_csv)


def reinvent(
    selection: str = "sele",
    prior: str = "all",
    num_smiles: Any = 10000,
    temperature: Any = 1,
) -> Optional[str]:
    """Start REINVENT4 Mol2Mol sampling for exactly one ligand residue."""

    try:
        return _reinvent_dispatch(selection, prior, num_smiles, temperature)
    except Exception as exc:
        print(f"reinvent: Error: {exc}")
        return None


# Visualization helpers
one_letter = {
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "GLU": "E",
    "GLN": "Q",
    "ASP": "D",
    "ASN": "N",
    "HIS": "H",
    "TRP": "W",
    "PHE": "F",
    "TYR": "Y",
    "ARG": "R",
    "LYS": "K",
    "SER": "S",
    "THR": "T",
    "MET": "M",
    "ALA": "A",
    "GLY": "G",
    "PRO": "P",
    "CYS": "C",
}

cmd.set("label_position", (3, 1, 1))
cmd.set("label_size", 25)
cmd.alias(
    "tag",
    'label n. CA and sele, "%s%s" % (one_letter[resn], resi)',
)
cmd.alias("b", "label n. CA and sele, b")
cmd.alias("untag", "hide label, sele")

_RESIDUE_STICK_ALIASES = {
    "k": "lys",
    "c": "cys",
    "s": "ser",
    "h": "his",
    "t": "thr",
    "r": "arg",
    "d": "asp",
    "e": "glu",
    "n": "asn",
    "q": "gln",
    "y": "tyr",
    "f": "phe",
    "w": "trp",
    "g": "gly",
    "a": "ala",
    "m": "met",
    "p": "pro",
    "i": "ile",
    "l": "leu",
    "v": "val",
}
for alias_name, residue_name in _RESIDUE_STICK_ALIASES.items():
    cmd.alias(
        alias_name,
        "hide sticks, polymer.protein; "
        f"show sticks, enabled and resn {residue_name}",
    )

cmd.alias(
    "me",
    "hide sticks, polymer.protein; "
    "show sticks, enabled and resn val; "
    "show sticks, enabled and resn leu; "
    "show sticks, enabled and resn ile; "
    "show sticks, enabled and resn met",
)
cmd.alias("z", "hide sticks, polymer.protein")
cmd.alias("addh", "h_add sele")
cmd.alias("remh", "remove sele extend 1 & hydro & not nbr. (don.|acc.)")
cmd.alias("del", "as nb_spheres, solvent; hide nb_spheres")
cmd.alias("sol", "as nb_spheres, solvent")
cmd.alias("ph", "hide (hydro) and enabled and (elem C extend 1)")


def white() -> None:
    """Color only enabled protein polymers white by atom."""

    util.cba(144, "enabled and polymer.protein", _self=cmd)


cmd.extend("white", white)


def cc() -> None:
    """Color enabled protein chains by chain identifier."""

    util.color_chains("(enabled and name CA)", _self=cmd)


cmd.extend("cc", cc)


def dssp() -> None:
    """Color enabled proteins by secondary structure."""

    cmd.color("gray", "enabled and ss l")
    cmd.color("yellow", "enabled and ss s")
    cmd.color("red", "enabled and ss h")


cmd.extend("dssp", dssp)


def plddt() -> None:
    """Color enabled AlphaFold structures by pLDDT stored in B-factors."""

    cmd.color("0x0053D6", "enabled and b < 100")
    cmd.color("0x65CBF3", "enabled and b < 90")
    cmd.color("0xFFDB13", "enabled and b < 70")
    cmd.color("0xFF7D45", "enabled and b < 50")
    util.cnc("enabled", _self=cmd)
    cmd.color("white", "enabled and sc. and elem C")


cmd.extend("plddt", plddt)


def _distance_text(distance: Any, command_name: str) -> Optional[str]:
    try:
        value = float(distance)
    except (TypeError, ValueError):
        print(f"{command_name}: Error: Distance must be a number.")
        return None
    if value < 0:
        print(f"{command_name}: Error: Distance must not be negative.")
        return None
    return f"{value:g}"


def see(distance: Any = 8) -> None:
    """Zoom to enabled ligand sites within a distance in angstroms."""

    distance_value = _distance_text(distance, "see")
    if distance_value is not None:
        cmd.zoom(f"byres organic around {distance_value} and enabled")


cmd.extend("see", see)


def line(distance: Any = 4) -> None:
    """Show lines for enabled residues near organic ligands."""

    distance_value = _distance_text(distance, "line")
    if distance_value is not None:
        cmd.show("lines", f"byres organic around {distance_value} and enabled")


cmd.extend("line", line)


def grab(sele: str = "organic", distance: Any = 4) -> None:
    """Select enabled residues near a ligand selection."""

    selection = str(sele).strip() or "organic"
    distance_value = _distance_text(distance, "grab")
    if distance_value is None:
        return
    cmd.select(
        "sele",
        f"byres ({selection}) around {distance_value} and enabled",
    )
    cmd.enable("sele")


cmd.extend("grab", grab)


def hb() -> None:
    """Show state-aware polar contacts from visible HETATM ligands."""

    ligand_selection = "visible and hetatm and not solvent"
    receptor_selection = "enabled and polymer.protein"
    water_selection = "enabled and solvent"

    if cmd.count_atoms(ligand_selection) == 0:
        print("hb: Error: No visible non-water HETATM ligand was found.")
        return

    cmd.delete("contact")
    cmd.delete("ligand")
    cmd.delete("hb_receptor")
    cmd.delete("hb_water")

    created_objects = []
    if cmd.count_atoms(receptor_selection) > 0:
        cmd.distance(
            "hb_receptor",
            ligand_selection,
            receptor_selection,
            mode=2,
            label=0,
            reset=1,
            state=0,
            quiet=1,
        )
        cmd.enable("hb_receptor")
        created_objects.append("hb_receptor")

    if cmd.count_atoms(water_selection) > 0:
        cmd.distance(
            "hb_water",
            ligand_selection,
            water_selection,
            mode=2,
            label=0,
            reset=1,
            state=0,
            quiet=1,
        )
        cmd.enable("hb_water")
        created_objects.append("hb_water")

    if not created_objects:
        print("hb: Error: No enabled protein receptor or water atoms were found.")
        return

    print(
        "hb: Created state-aware polar contacts: "
        + ", ".join(created_objects)
        + "."
    )


cmd.extend("hb", hb)


def dock(ligand: str = "enabled and organic") -> None:
    """Show polar contacts between a ligand selection and other objects."""

    ligand_selection = str(ligand).strip() or "enabled and organic"
    cmd.delete("polar_conts")
    cmd.dist(
        "polar_conts",
        f"({ligand_selection})",
        f"not ({ligand_selection})",
        quiet=1,
        mode=2,
        label=0,
        reset=1,
    )
    cmd.enable("polar_conts")


cmd.extend("dock", dock)


def at(selection: str = "sele") -> None:
    """Show state-aware protein and water contacts for a selection."""

    target_selection = str(selection).strip() or "sele"
    if cmd.count_atoms(target_selection) == 0:
        print(f"at: Error: No atoms were found in selection '{target_selection}'.")
        return

    receptor_selection = (
        f"(enabled and polymer.protein) and not ({target_selection})"
    )
    water_selection = f"(enabled and solvent) and not ({target_selection})"

    cmd.delete("contact")
    cmd.delete("at_receptor")
    cmd.delete("at_water")

    created_objects = []
    if cmd.count_atoms(receptor_selection) > 0:
        cmd.distance(
            "at_receptor",
            target_selection,
            receptor_selection,
            quiet=1,
            mode=2,
            label=0,
            reset=1,
            state=0,
        )
        cmd.enable("at_receptor")
        created_objects.append("at_receptor")

    if cmd.count_atoms(water_selection) > 0:
        cmd.distance(
            "at_water",
            target_selection,
            water_selection,
            quiet=1,
            mode=2,
            label=0,
            reset=1,
            state=0,
        )
        cmd.enable("at_water")
        created_objects.append("at_water")

    if not created_objects:
        print("at: Error: No enabled protein receptor or water atoms were found.")
        return

    print(
        "at: Created state-aware polar contacts: "
        + ", ".join(created_objects)
        + "."
    )


cmd.extend("at", at)


def byres(ligand: str = "enabled and organic", distance: Any = 10) -> None:
    """Create ligand and nearby-residue objects for interaction analysis."""

    distance_value = _distance_text(distance, "byres")
    if distance_value is None:
        return
    ligand_selection = str(ligand).strip() or "enabled and organic"
    if cmd.count_atoms(ligand_selection) == 0:
        print(f"byres: Error: No atoms were found in selection '{ligand_selection}'.")
        return

    object_suffix = distance_value.replace("-", "m").replace(".", "_")
    object_name = f"byres{object_suffix}"
    nearby_selection = (
        f"byres ((enabled and not ({ligand_selection})) within "
        f"{distance_value} of ({ligand_selection}))"
    )
    source_selection = f"({ligand_selection}) or ({nearby_selection})"
    source_objects = list(cmd.get_object_list(source_selection))

    cmd.create("Lig", ligand_selection)
    cmd.create(object_name, nearby_selection)

    for source_object in source_objects:
        if source_object not in {"Lig", object_name}:
            cmd.disable(source_object)

    cmd.hide("everything", "Lig")
    cmd.hide("everything", object_name)
    cmd.show("sticks", "Lig")
    line_selection = (
        f"byres (({object_name} and polymer.protein) within 4 of Lig)"
    )
    cmd.show("lines", line_selection)

    cmd.enable("Lig")
    cmd.enable(object_name)
    cmd.zoom(f"(Lig) or ({line_selection})", buffer=2.0)
    hb()
    cmd.do("sol")


cmd.extend("byres", byres)


def well(distance: Any = 5) -> None:
    """Show the surface of enabled residues near organic ligands."""

    distance_value = _distance_text(distance, "well")
    if distance_value is not None:
        cmd.show("surface", f"byres organic expand {distance_value} and enabled")


cmd.extend("well", well)


def cav() -> None:
    """Enable PyMOL surface cavity mode."""

    cmd.set("surface_cavity_mode", 2)


cmd.extend("cav", cav)


def wall() -> None:
    """Disable PyMOL surface cavity mode."""

    cmd.set("surface_cavity_mode", 0)


cmd.extend("wall", wall)


def pl() -> None:
    """Apply the ligand-cartoon preset to enabled objects."""

    preset.ligand_cartoon("enabled", _self=cmd)
    cmd.show("spheres", "solvent")


cmd.extend("pl", pl)


_DOCKING_SCORE_NAME = re.compile(
    r"(?:score|affinity|energy|confidence|binding|rmsd)",
    flags=re.IGNORECASE,
)
_NON_SCORE_METRIC_NAME = re.compile(
    r"(?:rank|pose|runtime|time)",
    flags=re.IGNORECASE,
)
_MAX_SCORE_SDF_BYTES = 100 * 1024 * 1024


def _object_property_names(object_name: str, state: int) -> list[str]:
    """Return normalized object-property names for one PyMOL state."""

    raw_properties = cmd.get_property_list(object_name, state=state) or []
    if isinstance(raw_properties, dict):
        raw_properties = raw_properties.keys()

    names = []
    for item in raw_properties:
        if isinstance(item, str):
            property_name = item
        elif isinstance(item, (list, tuple)) and item:
            property_name = str(item[0])
        else:
            continue
        if property_name and property_name not in names:
            names.append(property_name)
    return names


def _sdf_record_properties(record: str) -> dict[str, str]:
    """Parse SD tags from one SDF record without requiring RDKit."""

    properties = {}
    lines = record.splitlines()
    index = 0
    while index < len(lines):
        match = re.match(r"^>\s*<([^>]+)>", lines[index])
        if not match:
            index += 1
            continue
        property_name = match.group(1).strip()
        index += 1
        values = []
        while index < len(lines) and lines[index].strip():
            values.append(lines[index])
            index += 1
        properties[property_name] = "\n".join(values).strip()
    return properties


def _sdf_state_properties(sdf_text: str, state: int) -> dict[str, str]:
    """Return SD tags for a one-based record/state index."""

    records = [
        record
        for record in re.split(r"^\$\$\$\$\s*$", sdf_text, flags=re.MULTILINE)
        if record.strip()
    ]
    if state < 1 or state > len(records):
        return {}
    return _sdf_record_properties(records[state - 1])


def _read_sdf_state_properties(path: Path, state: int) -> dict[str, str]:
    """Read one state's SD tags from a size-limited local SDF file."""

    resolved_path = path.expanduser().resolve()
    if not resolved_path.is_file():
        return {}
    if resolved_path.stat().st_size > _MAX_SCORE_SDF_BYTES:
        raise ValueError("the SDF file exceeds the 100 MiB size limit")
    sdf_text = resolved_path.read_text(encoding="utf-8", errors="replace")
    return _sdf_state_properties(sdf_text, state)


def _score_source_candidates(object_name: str, source_file: str) -> list[Path]:
    """Return explicit or narrowly inferred SDF source paths."""

    if source_file:
        return [Path(source_file)]
    filename = f"{object_name}.sdf"
    candidates = [Path.cwd() / filename, Path.home() / "Downloads" / filename]
    return list(dict.fromkeys(candidates))


def _state_properties(
    object_name: str,
    state: int,
    source_file: str,
) -> dict[str, Any]:
    """Read state properties through modern and legacy PyMOL fallbacks."""

    get_property_list = getattr(cmd, "get_property_list", None)
    get_property = getattr(cmd, "get_property", None)
    if callable(get_property_list) and callable(get_property):
        property_names = _object_property_names(object_name, state)
        properties = {
            name: get_property(name, object_name, state=state)
            for name in property_names
        }
        if properties:
            return properties

    get_sdfstr = getattr(cmd, "get_sdfstr", None)
    if callable(get_sdfstr):
        try:
            exported_sdf = get_sdfstr(f"({object_name})", state=state, quiet=1)
        except Exception:
            exported_sdf = None
        if exported_sdf:
            properties = _sdf_state_properties(str(exported_sdf), 1)
            if properties:
                return properties

    for candidate in _score_source_candidates(object_name, source_file):
        properties = _read_sdf_state_properties(candidate, state)
        if properties:
            return properties
    return {}


def _numeric_property_value(value: Any) -> Optional[Any]:
    """Return a numeric property value, accepting numeric SDF strings."""

    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    try:
        numeric_value = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return int(numeric_value) if numeric_value.is_integer() else numeric_value


def _deduplicate_score_aliases(scores: dict[str, Any]) -> dict[str, Any]:
    """Remove equal ``Affinity``/``Score_affinity`` aliases."""

    deduplicated = {}
    aliases: dict[str, str] = {}
    for property_name, value in scores.items():
        alias = re.sub(r"(^|_)score_", r"\1", property_name.casefold())
        existing_name = aliases.get(alias)
        if existing_name is None or deduplicated[existing_name] != value:
            deduplicated[property_name] = value
            aliases[alias] = property_name
            continue

        existing_is_score = "_score_" in existing_name.casefold()
        current_is_score = "_score_" in property_name.casefold()
        if current_is_score and not existing_is_score:
            del deduplicated[existing_name]
            deduplicated[property_name] = value
            aliases[alias] = property_name
    return deduplicated


def printscore(
    selection: str = "enabled",
    score_name: str = "",
    source_file: str = "",
) -> dict[str, dict[str, Any]]:
    """Print docking scores for the current state of selected objects."""

    target_selection = str(selection).strip() or "enabled"
    requested_name = str(score_name).strip()
    if requested_name.casefold() in {"auto", "*"}:
        requested_name = ""
    source_path = str(source_file).strip()
    try:
        object_names = list(cmd.get_object_list(f"({target_selection})"))
    except Exception as exc:
        print(f"printscore: Error: Invalid selection '{target_selection}': {exc}")
        return {}

    if not object_names:
        print(f"printscore: Error: No molecular objects match '{target_selection}'.")
        return {}

    current_state = int(cmd.get_state())
    results: dict[str, dict[str, Any]] = {}
    inspected_objects = []

    for object_name in object_names:
        state_count = int(cmd.count_states(object_name))
        if state_count < current_state:
            continue
        inspected_objects.append(object_name)

        try:
            properties = _state_properties(object_name, current_state, source_path)
        except Exception as exc:
            print(
                f"printscore: Warning: Could not read properties from "
                f"'{object_name}' state {current_state}: {exc}"
            )
            continue
        property_names = list(properties)

        if requested_name:
            matching_names = [
                name
                for name in property_names
                if name.casefold() == requested_name.casefold()
            ]
        else:
            matching_names = [
                name
                for name in property_names
                if _DOCKING_SCORE_NAME.search(name)
                and not _NON_SCORE_METRIC_NAME.search(name)
            ]

        object_scores = {}
        for property_name in matching_names:
            raw_value = properties[property_name]
            value = _numeric_property_value(raw_value)
            if value is not None:
                object_scores[property_name] = value
        if not requested_name:
            object_scores = _deduplicate_score_aliases(object_scores)

        if object_scores:
            results[object_name] = object_scores
            for property_name, value in object_scores.items():
                print(
                    f"printscore: object={object_name}, state={current_state}, "
                    f"{property_name}={value}"
                )

    if not results:
        if not inspected_objects:
            print(
                f"printscore: Error: State {current_state} is not available in "
                "the matching objects."
            )
        elif requested_name:
            print(
                f"printscore: Error: Numeric property '{requested_name}' was not "
                f"found in state {current_state}."
            )
        else:
            print(
                f"printscore: Error: No numeric docking score properties were "
                f"found in state {current_state}. If this PyMOL version does "
                "not preserve SDF properties, provide the original SDF as the "
                "third argument."
            )
    return results


cmd.extend("printscore", printscore)


def _print_and_copy_fastastr(command_name: str, selection: str) -> str:
    sequence = cmd.get_fastastr(selection)
    print(sequence)
    try:
        _copy_to_clipboard(sequence)
        print(f"{command_name}: Copied the FASTA sequence to the clipboard.")
    except SmilesError as exc:
        print(f"{command_name}: Warning: Could not copy to the clipboard: {exc}")
    return sequence


def fasta() -> str:
    """Print, copy, and return the FASTA sequence of enabled objects."""

    return _print_and_copy_fastastr("fasta", "enabled")


cmd.extend("fasta", fasta)


def seq() -> str:
    """Print, copy, and return the FASTA sequence of ``sele``."""

    return _print_and_copy_fastastr("seq", "sele")


cmd.extend("seq", seq)


# Register these after all single-letter aliases so the full command names
# always have the final entries in PyMOL's command table.
cmd.extend("reinvent_priors", reinvent_priors)
cmd.extend("reinvent", reinvent)
print("chemistry: Loaded all commands, including reinvent and reinvent_priors.")
