from matplotlib import pyplot as plt
import numpy as np
import sys
import subprocess
# TORUN: python save_colormap.py acton = colorname from the library
# Try importing colormaps, install if missing
try:
    import colormaps as cmaps
except ImportError:
    print("The 'colormaps' package is not installed. Installing now...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "colormaps"])  # colorcet contains these maps
    import colormaps as cmaps  # fallback to colorcet

def save_colormap(cmap_name, filename='colormap.txt', n=256):
    try:
        # Retrieve colormap from colormaps or colorcet
        cmap_func = getattr(cmaps, cmap_name)
        cmap = plt.get_cmap(cmap_func)
    except AttributeError:
        print(f"Colormap '{cmap_name}' not found in colormaps module.")
        sys.exit(1)

    # Sample the colormap
    colors = cmap(np.linspace(0, 1, n))[:, :3]  # Drop alpha channel

    # Save to file
    with open(filename, 'w') as f:
        f.write('[\\n')
        for row in colors:
            f.write('  {:.6f} {:.6f} {:.6f};\\n'.format(*row))
        f.write(']\\n')

    # Plot the colormap
    plt.imshow([colors], aspect='auto')
    plt.axis('off')
    plt.title(cmap_name)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python save_cmap.py <colormap_name> [output_filename]")
        sys.exit(1)

    cmap_name = sys.argv[1]
    output_filename = sys.argv[2] if len(sys.argv) > 2 else f"{cmap_name}.txt"

    save_colormap(cmap_name, output_filename)

