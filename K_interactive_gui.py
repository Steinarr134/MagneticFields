import matplotlib.pyplot as plt

from K_interactive_simulations import *
from matplotlib.widgets import Button, Slider, RadioButtons
from matplotlib import gridspec

plt.style.use('dark_background')
# plt.rcParams['text.usetex'] = True

# Create the figure and the line that we will manipulate
fig1, ax1 = plt.subplots()
fig1.subplots_adjust(left=0.02, right=0.98, bottom=0.1, top=0.98, wspace=0.4)
ax1.remove()
# gs1 = fig1.add_gridspec(2, 3)
# Axes["Current"] = fig1.add_subplot(gs1[0, 1:])
# Axes["Ellipse"] = fig1.add_subplot(gs1[:, 0])
# Axes["Sequence"] = fig1.add_subplot(gs1[1, 1:])

gs1 = fig1.add_gridspec(3, 7)
Axes["Current"] = fig1.add_subplot(gs1[0, 1:3])
Axes["Ellipse"] = fig1.add_subplot(gs1[:2, 0])
Axes["Sequence"] = fig1.add_subplot(gs1[1, 1:3])

# place sliders in their own sub-gridspec, use 8 columns and leave the first onw empty
# to make room for the labels
sliderspec = gridspec.GridSpecFromSubplotSpec(8, len(params1), gs1[2, :2])

for i, param in enumerate(params1):
    params1[param]['ax'] = plt.subplot(sliderspec[1:, i])
    params1[param]['slider'] = Slider(
        ax=params1[param]['ax'],
        label=params1[param]['dispname'],
        valmin=params1[param]['min'],
        valmax=params1[param]['max'],
        valinit=params1[param]['init'],
        orientation="vertical"
    )
    params1[param]['slider'].on_changed(lambda val, param=param: update(param, val))

Axes["text"] = [fig1.add_subplot(gs1[2, 2:4])]
Axes["text"][0].xaxis.set_visible(False)
Axes["text"][0].yaxis.set_visible(False)
Axes["text"].append(Axes["text"][0].text(0.1, 0.3, "hellllllloooooo"))

ResetButton = Button(fig1.add_axes([0.95, 0.01, 0.05, 0.02]), "reset")
def reset_params(_):
    for param in params1:
        params1[param]['slider'].reset()
ResetButton.on_clicked(reset_params)
# # The function to be called anytime a slider's value changes
# def update(val):
#     line.set_ydata(f(t, amp_slider.val, freq_slider.val))
#     fig.canvas.draw_idle()


# # register the update function with each slider
# freq_slider.on_changed(update)
# amp_slider.on_changed(update)

# create second figure that shows B1/B2 as a function of thing to be chosen

# fig2, ax2 = plt.subplots()
# ax2.remove()
# gs1 = fig2.add_gridspec(1,1)
# Axes["chosen"] = fig2.add_subplot(gs1[0, 0])
# fig2.subplots_adjust(left=0.3)
# Axes["Radio"] = [fig2.add_axes([0.02, 0.5, 0.1, 0.4])]

Axes["chosen"] = fig1.add_subplot(gs1[:, 4:])
Axes["Radio"] = [fig1.add_subplot(gs1[0:2, 3])]

radio = RadioButtons(Axes["Radio"][0], [params1[p]['dispname'] for p in params1])
radio.on_clicked(plot_as_a_function_of)
Axes["Radio"].append(radio)
# Create selection


update("D", 0)
plt.show()
