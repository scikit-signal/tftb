from tftb.processing.linear import ShortTimeFourierTransform
from tftb.generators import fmconst
import numpy as np
sig = np.r_[fmconst(128, 0.2)[0], fmconst(128, 0.4)[0]]
tfr = ShortTimeFourierTransform(sig)
tfr.run()
tfr.plot()
