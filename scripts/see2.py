#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


def main():
    nrows = 3
    fig, axes = plt.subplots(nrows, 2)

    for row in axes:
        x = np.random.normal(0, 1, 100).cumsum()
        y = np.random.normal(0, 0.5, 100).cumsum()
        plot(row, x, y)

    plt.show()


def plot(axrow, x, y):
    axrow[0].plot(x, color='red')
    axrow[1].plot(y, color='green')


main()
