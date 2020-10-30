import matplotlib.pyplot as plt
import numpy as np

from fompy.functions import fd1

if __name__ == '__main__':
    xs = np.linspace(-4, 1, 1000)

    k = fd1(-10) / np.exp(-10)

    plt.plot(xs, fd1(xs), label='F(1/2;x)')
    plt.plot(xs, np.exp(xs), label='exp(x)')
    plt.plot(xs, np.exp(xs) * k, label='exp(x)*k')
    plt.grid()
    plt.legend()
    plt.show()
