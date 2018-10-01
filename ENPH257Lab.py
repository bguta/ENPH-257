import serial
import numpy as np
import matplotlib.pyplot as plt

port = "COM6"
b_animate = True

L = 2.54 * 12  # cm
T_amb = 20
Sensors = 5

T = np.zeros((1, 5))
x = np.ones((1, 5))

# distances for each sensor
x[0] = 1.3 / 100
x[1] = 8.3 / 100
x[2] = 15.3 / 100
x[3] = 22.2 / 100
x[4] = 29.25 / 100


def show():


if b_animate:
    # initialize plot
    fig = plt.figure()
    ax = plt.axes([0, 10, 20, 80])
    ax.grid(True)
    plt.xlabel('time $m$', fontsize=10)  # X axis label
    plt.ylabel('Temperature in $^o C$', fontsize=10)  # Y axis label
    plt.title('Tempurature in an Aluminum Rod\nheated by $10W$ in a 25$^o C$ Environment')
    plt.legend()
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    temp_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    line, = ax.plot([], [], lw=2, color='xkcd:salmon')

    def init():
        '''
        init function required for animator
        '''
        line.set_data(x, u_1)
        time_text.set_text('')
        temp_text.set_text('')

        return line, time_text, temp_text

    def animate(i):
        global et, dt, done, line, T

        T =  # function for temp

        # increment time
        et += dt

        # Update changing info for the animator
        line.set_data(x, u_1)
        time_text.set_text('time = {0:.4f}s'.format(et))
        temp_text.set_text('Average Tempurature = {0:.2f}$^o C$'.format(sum(u_1) / len(u_1)))

        # Return a tuple of the updated fields for the animator
        return line, time_text, temp_text

    #Animate and plot
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=Nt, interval=dt / 50, blit=True)
