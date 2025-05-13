import py3Dmol
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def show_results(cif_file):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(open(cif_file).read(), 'cif')

    view.setStyle({}, {'cartoon': {}})

    bins = [('90-100', '#0053D6'), ('70-90',  '#65CBF3'), ('50-70',  '#FFDB13'), ('0-50',   '#FF7D45'), ]
    for rng, col in bins:
        view.setStyle({'b': rng}, {'cartoon': {'color': col}})

    view.setStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})

    view.setBackgroundColor('white')
    view.zoomTo()
    return view


def create_pae_plot(
    data, model_name="", cmap="Greens_r", vmin=0, vmax=31.75
):
    """Create a PAE plot from data matrix.

    Args:
        data (dict): Dictionary containing PAE data
        model_name (str): Name of the model to display in the title
        cmap (str): Colormap for the plot
        vmin (float): Minimum value for colormap
        vmax (float): Maximum value for colormap
    """
    fig, ax = plt.subplots(figsize=(3.6, 4.2))

    if model_name:
        ax.set_title(model_name)

    ax.set_xlabel("Scored Residue")
    ax.set_ylabel("Aligned Residue")

    mappable = ax.imshow(
        data["pae"],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    ax.set_xlim(0, data["pae"].shape[1])
    ax.set_ylim(data["pae"].shape[0], 0)
    ax.set_aspect("equal")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.5)
    cbar = fig.colorbar(
        mappable,
        cax=cax,
        ax=ax,
        orientation="horizontal",
        pad=0.2,
        label="Expected Position Error (Ångströms)"
    )
