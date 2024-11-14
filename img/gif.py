#%%
import os
import imageio.v2 as imageio  # 変更点はこちら

# スクリプトが置かれているディレクトリを取得
script_dir = os.path.dirname(os.path.abspath(__file__))

# .pngファイルを取得し、ソート
image_files = sorted([os.path.join(script_dir, f) for f in os.listdir(script_dir) if f.endswith('.png')])

# GIFとして保存（ループ設定あり）
gif_filename = os.path.basename(script_dir) + '.gif'
with imageio.get_writer(gif_filename, mode='I', duration=1000, loop=0) as writer:
    for filename in image_files:
        image = imageio.imread(filename)
        writer.append_data(image)

print(f"ループするGIFが作成されました: {gif_filename}")
