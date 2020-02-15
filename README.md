# MA5232 Project 1 HUANG LINHANG
The repository is dedicated to the project for the NUS module MA5232 coordinated by Professor Bao Weizhu. The codes might not be well-commented as they only serves for the purpose of answering assigned questions. The author is not responsible for the potential consequences of implementing the codes. Note that the author was writing the algorithms without enough knowledge of relavant, well-studied libraries in Julia. Therefore, it is highly likely that the codes are not optimised very well. Thanks for understanding.

1. The codes of imaginary-time algorithms for 1D and 2D situation can be referred from ***imag_time.jl***
   - Note that the author has used default linear solver in Julia, which is not optimal in time-complexity for 2D situation. Hence, some modification is suggested for high-resolution 2D implementation.
   - Only the cases of squared range in 2D is included.

2. The codes of time-splitting spectral method for 1D and 2D situation can be referred from ***tssp.jl***
   - The author only has designed the codes up to second-order precision.
   - Just like the imaginary-time method, the author only has considered the cases of squared range in 2D.

3. The codes of implementing 1D tssp with special initial condition can be found in ***dynamics.jl***
