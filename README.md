# Hilbert transform

fftw3 is required.

The function's output is same with scipy which is a digital process package of python.

### hilbert_r2c()

* `input`: double array
* `output`: fftw_complex array

### hilbert_r2r()

* `input`: double array
* `output`: double array

How to compare result with scipy:
```python
import scipy.signal.signaltools as sigtool
import numpy as np

N = 16
x = np.linspace(1, N, N)

print(sigtool.hilbert(x))   # same result with hilbert_r2c()
print(np.abs(sigtool.hilbert(x))) # same result with hilbert_r2r()
```
