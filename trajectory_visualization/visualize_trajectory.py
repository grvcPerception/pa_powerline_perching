# -*- coding: utf-8 -*-
"""
Drone animation adapted from:
https://github.com/bobzwik/Quadcopter_SimCon
author: John Bass
email: john.bobzwik@gmail.com
license: MIT

The rest of the code follows the MIT license of this project

Please feel free to use and modify this, but keep the above information. Thanks!
"""

import sys
import yaml
import argparse
import numpy as np
from numpy import pi
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.interpolate import interp1d
from pathlib import Path

def quat2Dcm(q):
    dcm = np.zeros([3,3])

    dcm[0,0] = q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
    dcm[0,1] = 2.0*(q[1]*q[2] - q[0]*q[3])
    dcm[0,2] = 2.0*(q[1]*q[3] + q[0]*q[2])
    dcm[1,0] = 2.0*(q[1]*q[2] + q[0]*q[3])
    dcm[1,1] = q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2
    dcm[1,2] = 2.0*(q[2]*q[3] - q[0]*q[1])
    dcm[2,0] = 2.0*(q[1]*q[3] - q[0]*q[2])
    dcm[2,1] = 2.0*(q[2]*q[3] + q[0]*q[1])
    dcm[2,2] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2

    return dcm

def sameAxisAnimation(t_all, pos_all, quat_all, vel_all, vel_norm, ifsave, vel_factor=1.0):

    # Init variables
    global quad_vel
    x = pos_all[:, 0]
    y = pos_all[:, 1]
    z = pos_all[:, 2]

    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax3d = fig.gca(projection='3d')
    projections = []
    quad_parts1, = ax.plot([], [], [], lw=2, color='red', zorder=5)
    quad_parts2, = ax.plot([], [], [], lw=2, color='blue', zorder=5)
    quad_traj, = ax.plot([], [], [], '--', lw=1, color='blue', zorder=5)

    # Setting the axes properties (all axes equal length)
    linesx =  np.array([ [lc_x[i]+lv_x[i]*sgm_length[i], lc_x[i]-lv_x[i]*sgm_length[i]] for i in range(len(lc_x))])
    linesy =  np.array([ [lc_y[i]+lv_y[i]*sgm_length[i], lc_y[i]-lv_y[i]*sgm_length[i]] for i in range(len(lc_y))])
    linesz =  np.array([ [lc_z[i]+lv_z[i]*sgm_length[i], lc_z[i]-lv_z[i]*sgm_length[i]] for i in range(len(lc_z))])
    extraEachSide = 0.5
    maxRange = 0.5*np.array([max(x.max(),linesx.max()) - min(x.min(),linesx.min()),
                             max(y.max(),linesy.max()) - min(y.min(),linesy.min()),
                             max(z.max(),linesz.max()) - min(z.min(),linesz.min())
                             ]).max() + extraEachSide
    mid_x = 0.5*(max(x.max(),linesx.max())+min(x.min(),linesx.min()))
    mid_y = 0.5*(max(y.max(),linesy.max())+min(y.min(),linesy.min()))
    mid_z = 0.5*(max(z.max(),linesz.max())+min(z.min(),linesz.min()))

    ax.set_xlim3d([mid_x-maxRange, mid_x+maxRange])
    ax.set_xlabel('X')
    ax.set_ylim3d([mid_y-maxRange, mid_y+maxRange])
    ax.set_ylabel('Y')
    ax.set_zlim3d([mid_z-maxRange, mid_z+maxRange])
    ax.set_zlabel('Z')

    # Set lines
    line1, = ax.plot([lc_x[0]-lv_x[0]*sgm_length[0], lc_x[0]+lv_x[0]*sgm_length[0]],
                     [lc_y[0]-lv_y[0]*sgm_length[0], lc_y[0]+lv_y[0]*sgm_length[0]], 
                     [lc_z[0]-lv_z[0]*sgm_length[0], lc_z[0]+lv_z[0]*sgm_length[0]], '-', lw=3, color='green', zorder=10)
    line2, = ax.plot([lc_x[1]-lv_x[1]*sgm_length[1], lc_x[1]+lv_x[1]*sgm_length[1]],
                     [lc_y[1]-lv_y[1]*sgm_length[1], lc_y[1]+lv_y[1]*sgm_length[1]], 
                     [lc_z[1]-lv_z[1]*sgm_length[1], lc_z[1]+lv_z[1]*sgm_length[1]], '-', lw=3, color='green', zorder=10)
    line3, = ax.plot([lc_x[2]-lv_x[2]*sgm_length[2], lc_x[2]+lv_x[2]*sgm_length[2]],
                     [lc_y[2]-lv_y[2]*sgm_length[2], lc_y[2]+lv_y[2]*sgm_length[2]], 
                     [lc_z[2]-lv_z[2]*sgm_length[2], lc_z[2]+lv_z[2]*sgm_length[2]], '-', lw=3, color='green', zorder=10)
    
    # Initial vel
    quad_vel = ax.quiver(pos_all[0][0], pos_all[0][1], pos_all[0][2], vel_all[0][0], vel_all[0][1], vel_all[0][2], length=vel_norm[0]*vel_factor)

    # Set time display
    titleTime = ax.text2D(0.05, 0.90, "", transform=ax.transAxes)

    def updateLines(i):

        global quad_vel

        time = t_all[i]
        pos = pos_all[i]
        vel = vel_all[i]
        vel_n = vel_norm[i]
        x = pos[0]
        y = pos[1]
        z = pos[2]

        x_from0 = pos_all[0:i, 0]
        y_from0 = pos_all[0:i, 1]
        z_from0 = pos_all[0:i, 2]

        # Draw new quad pose
        quat = quat_all[i]
        dzm = 0.01

        R = quat2Dcm(quat)
        motorPoints = np.array([[quad_lx[0], quad_ly[0], dzm], 
                                 [0, 0, 0],
                                 [quad_lx[3], quad_ly[3], dzm], 
                                 [quad_lx[1], quad_ly[1], dzm],
                                 [0, 0, 0],
                                 [quad_lx[2], quad_ly[2], dzm]])
        motorPoints = np.dot(R, np.transpose(motorPoints))
        motorPoints[0, :] += x
        motorPoints[1, :] += y
        motorPoints[2, :] += z

        quad_parts1.set_data(motorPoints[0, 0:3], motorPoints[1, 0:3])
        quad_parts1.set_3d_properties(motorPoints[2, 0:3])
        quad_parts2.set_data(motorPoints[0, 3:6], motorPoints[1, 3:6])
        quad_parts2.set_3d_properties(motorPoints[2, 3:6])
        quad_traj.set_data(x_from0, y_from0)
        quad_traj.set_3d_properties(z_from0)
        quad_vel.remove()
        quad_vel = ax.quiver(x, y, z, vel[0], vel[1], vel[2], length=vel_n*vel_factor)
        titleTime.set_text(u"Time = {:.2f} s\nVel = {:.2f} m/s".format(time,vel_n))

        return line3, line2, line1, quad_traj, quad_parts2, quad_parts1

    def ini_plot():

        quad_parts1.set_data(np.empty([1]), np.empty([1]))
        quad_parts1.set_3d_properties(np.empty([1]))
        quad_parts2.set_data(np.empty([1]), np.empty([1]))
        quad_parts2.set_3d_properties(np.empty([1]))
        quad_traj.set_data(np.empty([1]), np.empty([1]))
        quad_traj.set_3d_properties(np.empty([1]))

        return line3, line2, line1, quad_traj, quad_parts2, quad_parts1

    def on_click(event):
        azim, elev = ax3d.azim, ax3d.elev
        projections.append((azim, elev))
        print(azim, elev)

    # Create the Animation object
    line_ani = animation.FuncAnimation(fig, updateLines, init_func=ini_plot, frames=len(
        t_all), interval=20, blit=False)

    if(ifsave):
        line_ani.save((save_gif_name+'.gif'), dpi=80, writer='imagemagick', fps=25)

    # cid = fig.canvas.mpl_connect('button_release_event', on_click)
    try:
        plt.show()
    except:
        pass

    return line_ani


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Trajectory visualization script")
    parser.add_argument("-t", "--traj", help = "Trajectory file", required = True)
    parser.add_argument("-q", "--quad", help = "Quad file", required = True)
    parser.add_argument("-l", "--lines", help = "Lines file", required = True)
    parser.add_argument('--show_vel', help = "Show quad velocity overtime", dest='show_vel', action='store_true', default=False)
    parser.add_argument('--save', help = "Save animation as GIF", dest='save', action='store_true', default=False)
    argument = parser.parse_args()

    csv_file = argument.traj
    quad_path = argument.quad
    lines_path = argument.lines
    show_vel = argument.show_vel
    save_gif = argument.save
    save_gif_name = Path(csv_file).stem

    # Load YAML configurations
    quad_file = open(quad_path,'r')
    quad = yaml.load(quad_file, Loader=yaml.FullLoader)

    lines_file = open(lines_path,'r')
    lines = yaml.load(lines_file, Loader=yaml.FullLoader)
    
    lc_x = lines["lc_x"]
    lc_y = lines["lc_y"]
    lc_z = lines["lc_z"]
    lv_x = lines["lv_x"]
    lv_y = lines["lv_y"]
    lv_z = lines["lv_z"]
    rad_l = lines["rad_l"]
    sgm_length = lines["sgm_length"]

    quad_lx = [quad["lx_0"], quad["lx_1"], quad["lx_2"], quad["lx_3"]]
    quad_ly = [quad["ly_0"], quad["ly_1"], quad["ly_2"], quad["ly_3"]]

    # Load trajectory
    data = genfromtxt(csv_file, delimiter=',', skip_header=True)
    t_all = data[:, 0]
    pos_all = data[:, 1:4]
    quat_all = data[:, 4:8]
    vel_all = data[:, 8:11]

    # Interpolate every 20 ms
    fx = interp1d(t_all, pos_all[:,0], kind='cubic')
    fy = interp1d(t_all, pos_all[:,1], kind='cubic')
    fz = interp1d(t_all, pos_all[:,2], kind='cubic')
    fqw = interp1d(t_all, quat_all[:,0], kind='cubic')
    fqx = interp1d(t_all, quat_all[:,1], kind='cubic')
    fqy = interp1d(t_all, quat_all[:,2], kind='cubic')
    fqz = interp1d(t_all, quat_all[:,3], kind='cubic')
    fvx = interp1d(t_all, vel_all[:,0], kind='cubic')
    fvy = interp1d(t_all, vel_all[:,1], kind='cubic')
    fvz = interp1d(t_all, vel_all[:,2], kind='cubic')

    t_new = np.linspace(t_all[0], t_all[-1], num=int((t_all[-1]-t_all[0])*50+1), endpoint=True)
    pos_new = np.column_stack((fx(t_new), fy(t_new), fz(t_new)))
    quat_new = np.column_stack((fqw(t_new), fqx(t_new), fqy(t_new), fqz(t_new)))
    vel_new = np.column_stack((fvx(t_new), fvy(t_new), fvz(t_new)))
    vel_norm_max = np.amax(np.linalg.norm(vel_new,axis=1))
    vel_norm = np.linalg.norm(vel_new,axis=1)
    vel_new_normalized = vel_new/vel_norm_max

    # Plot
    ani = sameAxisAnimation(t_new, pos_new, quat_new, vel_new_normalized, vel_norm, save_gif, vel_factor=(0.25 if show_vel else 0))
