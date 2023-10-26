import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from PIL import Image
import io

input_file = '/rtmp/akumar/MERRA2_78vavrs/MERRA2_400.inst3_3d_asm_Nv.20230719.SUB.nc'

A = xr.open_dataset(input_file)
dataset = np.array([[A.isel(time=4)[var].isel(lev=lev_id).values for var in ('U', 'V', 'T', 'RH', 'OMEGA')] for lev_id in (2, 5, 7, 9)]).reshape(1, -1, 361, 576).squeeze()

images = []

for index in range(20):
    data = dataset[index, :, :].squeeze()
    plt.figure(figsize=(6, 4))
    plt.imshow(data, cmap='gist_rainbow')

    # Add a bold border
    for spine in plt.gca().spines.values():
        spine.set_linewidth(4)

    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()

    # Capture the current plot as an image
    img_buf = io.BytesIO()
    plt.savefig(img_buf, format='png')
    img_buf.seek(0)
    img = Image.open(img_buf)
    images.append(img)
    plt.close()

# Save the list of images as a GIF
images[0].save('output_var1.gif', save_all=True, append_images=images[1:], duration=300, loop=0)


