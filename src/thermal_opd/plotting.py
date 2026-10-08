"""Shared steady/transient temperature and refracted-ray field figures."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def save_temperature_opd_figure(path,r,phi,temperature,x,y,mask,opd,*,title,pv_m,continuous_pv_m=None):
    """Same 3 temperature planes + entry OPD, 61 levels and 180 dpi.

    Three temperatures share one scale within a snapshot. OPD is displayed
    relative to that aperture's minimum; the saved raw field is not altered.
    """
    plot_phi = np.append(phi,2*np.pi)
    tx,ty = r[:,None]*np.cos(plot_phi)*1e3,r[:,None]*np.sin(plot_phi)*1e3
    fig,axes = plt.subplots(2,2,figsize=(11,10),layout='constrained')
    low,high = float(np.min(temperature)),float(np.max(temperature))
    if high<=low:
        high = low+1e-12
    levels = np.linspace(low,high,61)
    for ax,field,name in zip(axes.flat,temperature,('Front: z=0','Mid-plane: z=L/2','Rear: z=L')):
        contour = ax.contourf(tx,ty,np.column_stack((field,field[:,0])),levels=levels,cmap='inferno')
        ax.set_title(name)
        fig.colorbar(contour,ax=ax,label='Temperature rise [K]',shrink=.85)
    optical = np.where(mask,(opd-np.min(opd[mask]))*1e9,np.nan)
    contour = axes.flat[-1].contourf(x*1e3,y*1e3,optical,levels=61,cmap='viridis')
    optical_title = f'Refracted-ray OPD, sampled PV={pv_m*1e9:.4f} nm'
    if continuous_pv_m is not None:
        optical_title += f'\nContinuous numerical PV={continuous_pv_m*1e9:.4f} nm'
    axes.flat[-1].set_title(optical_title)
    fig.colorbar(contour,ax=axes.flat[-1],label='OPD minus aperture minimum [nm]',shrink=.85)
    for ax in axes.flat:
        ax.set_aspect('equal');ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]')
    axes.flat[-1].set_xlabel('Entry x [mm]');axes.flat[-1].set_ylabel('Entry y [mm]')
    fig.suptitle(title)
    fig.savefig(Path(path),dpi=180)
    plt.close(fig)
