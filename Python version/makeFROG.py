from numpy import roll, size, fliplr, abs, ceil, outer, dot, conj
from numpy.fft import fft, fftshift
from matplotlib.pyplot import imshow

def makeFROG(electricField):

    N = size(electricField)

    print(electricField)

    # outer product form
    #TODO MEMORY ERROR
    electricFROG = outer(electricField, electricField)

    print(electricFROG)

    # row rotation
    for n in range(1, N):
        electricFROG[n, :] = roll(electricFROG[n, :], -n)

    # permute the columns to the right order
    electricFROG = fliplr(fftshift(electricFROG, 1))

    # FFT each column and put 0 frequency in the correct place:
    electricFROG = roll(fft(electricFROG, axis=1), int(ceil(N/2)-1))

    # generate FROG trace (= |field|^2)
    intensityFROG = abs(electricFROG)**2

    return (intensityFROG, electricFROG)
