import numpy as np
import matplotlib.pyplot as plt
from fdint import fdk

if __name__ == '__main__':
    xs = np.linspace(-4, 1, 1000)

    plt.plot(xs, fdk(0.5, xs), label='F(1/2;x)')
    plt.plot(xs, np.exp(xs), label='exp(x)')
    plt.grid()
    plt.legend()
    plt.show()
