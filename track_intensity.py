import tropycal.tracks as tracks
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"font.size": 17, "font.weight": "bold"})

#basin = tracks.TrackDataset(basin='north_atlantic',
#                            source='ibtracs', include_btk=False)

url = "https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2022-042723.txt"
basin = tracks.TrackDataset(basin='north_atlantic', atlantic_url=url)

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


plot_intensity(name='alicia', year=1983, start=2)
plt.savefig('../figures/intensity/Intensity_Alicia.jpeg')

plot_intensity(name='harvey', year=2017, start=2)
plt.savefig('../figures/intensity/Intensity_Harvey.jpeg')

plot_intensity(name='laura', year=2020, start=2)
plt.savefig('../figures/intensity/Intensity_Laura.jpeg')

plot_intensity(name='ike', year=2008, start=2)
plt.savefig('../figures/intensity/Intensity_Ike.jpeg')

plot_intensity(name='rita', year=2005, start=2)
plt.savefig('../figures/intensity/Intensity_Rita.jpeg')

plt.close()
