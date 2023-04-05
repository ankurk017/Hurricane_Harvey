# -*- coding: utf-8 -*-
"""track_intensity.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/193ZF_Q4km-LNwpAQVCZ3YrG8P6H0xJDg
"""


import tropycal.tracks as tracks
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"font.size": 17, "font.weight": "bold"})

basin = tracks.TrackDataset(basin='north_atlantic',
                            source='hurdat', include_btk=False)


def plot_intensity(name='michael', year=2018, start=1, skip=10):
    storm = basin.get_storm((name, year))

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    color = 'tab:blue'
    ax.plot(storm['date'][start:-skip], storm['vmax']
            [start:-skip], color=color, marker='d')
    ax.set_ylabel('Maximum Sustained Wind Speed (kt)', color=color)
    ax.set_xlabel('Date (YYYY-MM-DD)', color='k')
    ax.set_yticks(np.arange(0, 150, 30))
    ax.tick_params(axis='y', labelcolor=color)

    ax2 = ax.twinx()
    color = 'tab:red'
    ax2.plot(storm['date'][start:-skip], storm['mslp']
             [start:-skip], color=color, marker='d')
    ax2.set_ylabel('Minimum Sea Pressure Level (hPa)', color=color)
    ax.tick_params(axis='x', which='major', rotation=30, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    ax.grid()
    ax.set_title(f'{name.upper()} {year}')
    plt.tight_layout()


plot_intensity(name='ida', year=2021)
plt.savefig('../figures/IDA_intensity.jpeg')

plot_intensity(name='harvey', year=2017, start=20)
plt.savefig('../figures/Harvey_intensity.jpeg')

plot_intensity(name='nicholas', year=2021, skip=1)
plt.savefig('../figures/Nicholas_intensity.jpeg')

plt.close()
