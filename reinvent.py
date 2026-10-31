"""PyMOL Mol2Mol command extension. Load via: run /path/reinvent_portable.py"""
from __future__ import annotations

import csv
import datetime as dt
import json
import os
import shutil
from pathlib import Path
import re
import subprocess
import sys
import tempfile

# No machine-specific paths: determine the REINVENT environment at runtime.
# Optional overrides: REINVENT_PYTHON, REINVENT_PRIORS_DIR, REINVENT4_HOME.
def _script_path():
    # PyMOL `run` may overwrite __file__; co_filename retains the loaded filename.
    return Path(_script_path.__code__.co_filename).resolve()


def _python_candidates():
    candidates = []
    def add(path):
        if path:
            p = Path(path).expanduser()
            if p.is_file() and p not in candidates:
                candidates.append(p)
    add(os.environ.get('REINVENT_PYTHON'))
    add(sys.executable)
    for key in ('CONDA_PREFIX', 'VIRTUAL_ENV'):
        prefix = os.environ.get(key)
        if prefix:
            add(Path(prefix) / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python'))
    exe = shutil.which('reinvent')
    if exe:
        # Executable script's shebang usually identifies its environment.
        try:
            first = Path(exe).read_bytes().splitlines()[0].decode('utf-8')
            if first.startswith('#!') and not first.endswith('/env python'):
                add(first[2:].strip().split()[0])
        except (OSError, UnicodeError, IndexError):
            pass
        add(Path(exe).parent / ('python.exe' if os.name == 'nt' else 'python'))
    home = Path.home()
    for root in (home / 'miniforge3', home / 'miniconda3', home / 'anaconda3',
                 home / 'mambaforge', home / 'micromamba', home / '.conda'):
        for d in (root, *sorted((root / 'envs').glob('*'))):
            add(d / ('python.exe' if os.name == 'nt' else 'bin/python'))
    return candidates


def _runtime():
    # Test in subprocess to avoid importing a different torch into PyMOL.
    # Probe outside PyMOL's working directory: its output folder is named
    # 'reinvent' and can shadow an editable REINVENT4 installation.
    # Import the actual entry module rather than relying on reinvent.__file__,
    # which can be None for namespace packages.
    check = ('import pathlib,torch,rdkit,reinvent.Reinvent as entry; '
             'print(pathlib.Path(entry.__file__).resolve())')
    errors = []
    for python in _python_candidates():
        try:
            with tempfile.TemporaryDirectory(prefix='reinvent_probe_') as probe_dir:
                result = subprocess.run([str(python), '-c', check], cwd=probe_dir,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        text=True, timeout=30, check=False)
        except (OSError, subprocess.TimeoutExpired) as exc:
            errors.append(f'{python}: {exc}')
            continue
        if result.returncode != 0:
            errors.append(f'{python}: {result.stderr.strip()[-180:]}')
            continue
        origin = Path(result.stdout.strip().splitlines()[-1])
        if not origin.is_file() or origin.name != 'Reinvent.py':
            errors.append(f'{python}: invalid REINVENT4 module path: {origin}')
            continue
        # The console script is installed in the same Python environment.
        console = python.parent / ('reinvent.exe' if os.name == 'nt' else 'reinvent')
        if console.is_file():
            return {'python': str(python), 'executable': str(console),
                    'package': str(origin)}
        errors.append(f'{python}: REINVENT console command missing: {console}')
    raise RuntimeError('A working REINVENT4 Python environment was not found. '
                       'Install REINVENT4 with its dependencies, or set REINVENT_PYTHON '
                       'to the Python executable of that environment. Checked:\n' +
                       '\n'.join(errors[:12]))


def _prior_roots(runtime):
    home = Path.home()
    pkg = Path(runtime['package']).parent
    # runtime['package'] points to reinvent/Reinvent.py; pkg.parent is the repo.
    bases = [Path(s).expanduser() for key in ('REINVENT_PRIORS_DIR', 'REINVENT4_HOME')
             for s in os.environ.get(key, '').split(os.pathsep) if s]
    bases += [pkg.parent, pkg, _script_path().parent, Path.cwd(),
              Path.cwd().parent, Path(runtime['python']).parent.parent]
    # Conventional locations support non-editable installations.
    bases += [home / 'REINVENT4', home / 'reinvent4', home / 'Documents' / 'REINVENT4',
              home / 'Documents' / 'modeling' / 'REINVENT4',
              home / 'apps' / 'REINVENT4', home / 'projects' / 'REINVENT4']
    roots = []
    for base in bases:
        for root in (base, base / 'priors', base / 'models', base / 'checkpoints'):
            try:
                rp = root.resolve()
                if rp.is_dir() and rp not in roots:
                    roots.append(rp)
            except OSError:
                pass
    return roots



def _safe_filename(s):
    result = re.sub(r'[^A-Za-z0-9_-]+', '_', str(s)).strip('_')
    return result[:64] or 'UNKNOWN'


def _find_priors(runtime):
    files = {}
    for root in _prior_roots(runtime):
        # Look inside prior directories and common checkpoint locations, but
        # avoid crawling entire home directories or unrelated repositories.
        search_dirs = [root]
        search_dirs += [p for p in root.glob('*/priors') if p.is_dir()]
        search_dirs += [p for p in root.glob('*/models') if p.is_dir()]
        for directory in search_dirs:
            for p in sorted(directory.glob('*.prior')):
                if p.is_file() and 'mol2mol' in p.stem.lower():
                    files.setdefault(p.stem.lower(), p.resolve())
    return files


def reinvent_priors():
    """List Mol2Mol prior files discovered by this PyMOL extension."""
    runtime = _runtime()
    print(f'[REINVENT] Environment: {runtime["python"]}')
    print(f'[REINVENT] Package: {runtime["package"]}')
    print(f'[REINVENT] Loaded script: {_script_path()}')
    print('[REINVENT] Searched directories:')
    for root in _prior_roots(runtime):
        print('  ', root)
    priors = _find_priors(runtime)
    if priors:
        print('[REINVENT] Found Mol2Mol priors:')
        for name, path in priors.items():
            print(f'  {name}: {path}')
    else:
        print('[REINVENT] No Mol2Mol .prior files found.')
        print('Set a directory containing mol2mol*.prior in PyMOL using:')
        print('  import os; os.environ["REINVENT_PRIORS_DIR"] = "/absolute/path/to/priors"')
    return priors


def _select_priors(value, runtime):
    all_priors = _find_priors(runtime)
    if not all_priors:
        raise RuntimeError('No Mol2Mol priors found. Run reinvent_priors to inspect search paths. Set REINVENT_PRIORS_DIR within PyMOL to the folder containing Mol2Mol .prior models.')
    if not value or value.lower() in ('all', '*'):
        return list(all_priors.values())
    names = [s.strip() for s in value.split('+') if s.strip()]
    chosen = []
    for name in names:
        key = Path(name).stem.lower()
        matched = all_priors.get(key)
        if matched is None:
            candidates = [p for k, p in all_priors.items() if key in k]
            if len(candidates) == 1:
                matched = candidates[0]
        if matched is None:
            raise ValueError(f'Prior {name!r} not uniquely found. Available: {", ".join(all_priors)}')
        if matched not in chosen:
            chosen.append(matched)
    return chosen


def _detect_device():
    """Check a real tensor operation, not merely whether an MPS backend exists."""
    import torch
    if torch.cuda.is_available():
        try:
            x = torch.ones((4, 4), device='cuda')
            torch.cuda.synchronize()
            if (x @ x).device.type == 'cuda':
                return 'cuda'
        except (RuntimeError, AssertionError) as exc:
            print(f'[REINVENT] CUDA test failed: {exc}', flush=True)
    backend = getattr(torch.backends, 'mps', None)
    if backend is not None and backend.is_available():
        try:
            x = torch.ones((4, 4), device='mps')
            y = x @ x
            torch.mps.synchronize()
            if y.device.type == 'mps':
                return 'mps'
        except (RuntimeError, AssertionError, TypeError) as exc:
            print(f'[REINVENT] MPS test failed: {exc}', flush=True)
    return 'cpu'


def _reported_device(log_text):
    """Record what REINVENT reports, separately from the requested device."""
    matches = re.findall(r'Using (CPU|CUDA|MPS)\b[^\n]*', log_text, flags=re.IGNORECASE)
    return matches[-1].strip() if matches else 'not reported'


def _worker(manifest_file):
    from rdkit import Chem
    manifest_file = Path(manifest_file).resolve()
    manifest = json.loads(manifest_file.read_text(encoding='utf-8'))
    job_dir = manifest_file.parent
    manifest['device'] = _detect_device()
    manifest['device_report'] = {}
    manifest_file.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(f'[REINVENT] Auto-selected device: {manifest["device"]}', flush=True)
    raw_mol = Chem.MolFromMolFile(str(job_dir / 'selected_ligand.mol'), removeHs=False, sanitize=True)
    if raw_mol is None:
        raise RuntimeError('Ligand could not be parsed as a chemically valid MOL. Verify bond orders and charges in PyMOL.')
    mol = Chem.RemoveHs(raw_mol)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    if not smiles or Chem.MolFromSmiles(smiles) is None:
        raise RuntimeError('Could not obtain a valid SMILES from selected ligand.')
    (job_dir / 'input.smi').write_text(smiles + '\n', encoding='utf-8')
    (job_dir / 'input_smiles.txt').write_text(smiles + '\n', encoding='utf-8')

    command = Path(manifest['runtime']['executable'])
    if not command.is_file():
        raise RuntimeError(f'Missing REINVENT executable: {command}')

    merged_csv = Path(manifest['output_csv'])
    tmp_csv = merged_csv.with_suffix('.partial.csv')
    header = None
    total_rows = 0
    for prior_str in manifest['priors']:
        prior = Path(prior_str)
        tag = _safe_filename(prior.stem)
        target = job_dir / (tag + '.csv')
        config = job_dir / (tag + '.toml')
        import json as _json
        # JSON quoting produces valid TOML strings for these paths/values.
        config.write_text(
            'run_type = "sampling"\n'
            f'device = {_json.dumps(manifest["device"])}\n\n'
            '[parameters]\n'
            f'model_file = {_json.dumps(str(prior))}\n'
            f'smiles_file = {_json.dumps(str(job_dir / "input.smi"))}\n'
            f'output_file = {_json.dumps(str(target))}\n'
            f'num_smiles = {manifest["num_smiles"]}\n'
            f'temperature = {manifest["temperature"]}\n'
            'sample_strategy = "multinomial"\n'
            'unique_molecules = true\n'
            'randomize_smiles = true\n', encoding='utf-8'
        )
        print(f'[REINVENT] Starting {prior.name} on {manifest["device"]}', flush=True)
        run_env = os.environ.copy()
        if manifest['device'] == 'mps':
            run_env.setdefault('PYTORCH_ENABLE_MPS_FALLBACK', '1')
        result = subprocess.run([str(command), '--device', manifest['device'], str(config)], cwd=str(job_dir),
                                env=run_env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, errors='replace', check=False)
        (job_dir / (tag + '.log')).write_text(result.stdout or '', encoding='utf-8')
        reported = _reported_device(result.stdout or '')
        manifest['device_report'][tag] = {'requested': manifest['device'],
                                           'reinvent_log': reported,
                                           'exit_code': result.returncode}
        manifest_file.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        print(f'[REINVENT] Device: requested={manifest["device"]}, REINVENT log={reported}', flush=True)
        if manifest['device'] == 'mps' and reported.lower().startswith('using cpu'):
            print('[REINVENT] WARNING: REINVENT reports CPU despite explicit --device mps. '
                  'GPU computation is not confirmed; inspect the REINVENT version and device handler.', flush=True)
        if result.returncode != 0 or not target.is_file():
            print(f'[REINVENT] FAILED {tag}: exit={result.returncode}. See {tag}.log', flush=True)
            continue
        with target.open(newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            if not fields:
                print(f'[REINVENT] Empty CSV for {tag}', flush=True)
                continue
            if header is None:
                header = [c for c in fields if c != 'prior'] + ['prior']
                with tmp_csv.open('w', newline='', encoding='utf-8') as out:
                    csv.writer(out).writerow(header)
            count = 0
            with tmp_csv.open('a', newline='', encoding='utf-8') as out:
                writer = csv.DictWriter(out, fieldnames=header, extrasaction='ignore')
                for item in reader:
                    item['prior'] = prior.stem
                    writer.writerow(item)
                    count += 1
            total_rows += count
        print(f'[REINVENT] Finished {tag}: {count} rows', flush=True)
        # Maintain recoverable partial results if later priors fail.
        if tmp_csv.exists():
            import shutil
            shutil.copyfile(tmp_csv, merged_csv)
    if header is None:
        raise RuntimeError('All model runs failed. Inspect log files under ' + str(job_dir))
    print(f'[REINVENT] Complete: {merged_csv} ({total_rows} rows)', flush=True)


def _dispatch(selection='sele', prior='all', num_smiles=10000, temperature=1.0):
    from pymol import cmd
    number = int(num_smiles)
    temp = float(temperature)
    if number < 1 or not 0 < temp <= 10:
        raise ValueError('num_smiles must be positive; temperature must be > 0 and <= 10.')
    runtime = _runtime()
    chosen = _select_priors(str(prior), runtime)
    atoms = cmd.get_model(f'({selection})').atom
    if not atoms:
        raise ValueError(f'No atoms in selection: {selection}')
    # Avoid mistakenly converting a protein + ligand mixture or multiple ligands.
    residues = {(a.model, a.segi, a.chain, a.resi, a.resn) for a in atoms}
    if len(residues) != 1:
        raise ValueError('Select exactly one ligand residue (one object/chain/residue).')
    model, segi, chain, resi, resn = next(iter(residues))
    # Object names derived from PDB IDs are conventional, not guaranteed.
    match = re.search(r'(?<![a-z0-9])([0-9][a-z0-9]{3})(?![a-z0-9])', model.lower())
    pdb_id = match.group(1).upper() if match else _safe_filename(model)
    base_name = '_'.join((dt.datetime.now().strftime('%Y%m%d'),
                          _safe_filename(pdb_id), _safe_filename(resn)))
    output_dir = Path.cwd() / 'reinvent'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / (base_name + '.csv')
    job_dir = Path(tempfile.mkdtemp(prefix=base_name + '_', dir=output_dir))
    mol_path = job_dir / 'selected_ligand.mol'
    cmd.save(str(mol_path), selection, state=-1, format='mol')
    if not mol_path.exists() or mol_path.stat().st_size == 0:
        raise RuntimeError('PyMOL failed to save the selected molecule')
    manifest = job_dir / 'manifest.json'
    manifest.write_text(json.dumps({'device': 'auto', 'priors': [str(p) for p in chosen],
                                    'num_smiles': number, 'temperature': temp,
                                    'output_csv': str(output_csv), 'runtime': runtime}, indent=2), encoding='utf-8')
    worker_python = Path(runtime['python'])
    log_path = job_dir / 'job.log'
    with log_path.open('w', encoding='utf-8') as log:
        proc = subprocess.Popen([str(worker_python), str(_script_path()), '--worker', str(manifest)],
                                cwd=str(job_dir), stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    print(f'[REINVENT] Started PID={proc.pid}; priors={len(chosen)}; requested={number} per prior')
    print(f'[REINVENT] Output CSV: {output_csv}')
    print(f'[REINVENT] Logs: {log_path}')
    return str(output_csv)


def reinvent(selection='sele', prior='all', num_smiles=10000, temperature=1):
    """Run Mol2Mol sampling with tested CUDA > MPS > CPU selection.

    PyMOL: reinvent [selection [, prior [, num_smiles [, temperature]]]]
    """
    return _dispatch(selection, prior, num_smiles, temperature)


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != '--worker':
        raise SystemExit('Use this file inside PyMOL with: run pymol_reinvent.py')
    _worker(sys.argv[2])
else:
    from pymol import cmd
    cmd.extend('reinvent', reinvent)
    cmd.extend('reinvent_priors', reinvent_priors)
