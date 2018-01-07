from numpy import roll, size, transpose, fliplr, square
from numpy.fft import fft, ifft, ifftshift, fftshift

def makeFROG(electricField):

    N = size(electricField)
    electricFROG = electricField*transpose(electricField)

    # row rotation
    for n in range(1,N):
        electricFROG[n,:] = roll(electricFROG[n,:], -n)

    # permute the columns to the right order
    electricFROG = fliplr(fftshift(electricFROG, 2))

    # FFT each column and put 0 frequency in the correct place:
    electricFROG = roll(fft(electricFROG, axis=1), ceil(N/2)-1)

    # generate FROG trace (= |field|^2)
    intensityFROG = square(abs(electricFROG))

return (intensityFROG, electricFROG)
