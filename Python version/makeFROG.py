from numpy import roll, size, fliplr, abs, ceil
from numpy.fft import fft, ifft, ifftshift, fftshift

def makeFROG(electricField):

    N = size(electricField)

    # outer product form
    electricFROG = electricField@electricField.T

    # row rotation
    for n in range(1,N):
        electricFROG[n, :] = roll(electricFROG[n, :], -n)

    # permute the columns to the right order
    electricFROG = fliplr(fftshift(electricFROG, 2))

    # FFT each column and put 0 frequency in the correct place:
    electricFROG = roll(fft(electricFROG, axis=1), ceil(N/2)-1)

    # generate FROG trace (= |field|^2)
    intensityFROG = abs(electricFROG)**2

    return (intensityFROG, electricFROG)
