try:
    from CannonBallisticFunctions import make_dataset
except:
    from fallback_py_funcs.py_make_dataset import make_dataset

_make_dataset = make_dataset
def make_dataset(charges, length, max_height_above=256, max_height_below=256,
                 num_threads=16, verbose=True, max_steps=100000, max_length=600,
                 step=1, impossible_cutoff=50, max_delta_t_error=1,
                 amin=-30, amax=60, gravity=0.05, drag=0.99, num_iterations=5, num_elements=20,
                 check_impossible=True, lambertW = False):
    return _make_dataset(charges, length, max_height_above, max_height_below,
                         num_threads, verbose, max_steps,
                         max_length, step, impossible_cutoff,
                         max_delta_t_error,
                         amin, amax, gravity, drag, num_iterations, num_elements, check_impossible, lambertW)