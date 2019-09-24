import pyso3 as so3
import numpy as np

f_p = so3.create_parameter_dict(16, 8)
f_length = so3.f_size(f_p)
print(f_length)
f_before = np.random.normal(size=(f_length)) + 1j*np.random.normal(size=(f_length))

flmn = so3.forward(f_before, f_p)

f = so3.inverse(flmn, f_p)
flmn_new = so3.forward(f, f_p)

print(np.mean(np.abs(flmn_new-flmn)), np.max(np.abs(flmn_new-flmn)))
print(np.mean(np.abs(f_before-f)), np.max(np.abs(f_before-f)))


g_before = np.random.normal(size=(f_length)) + 1j*np.random.normal(size=(f_length))
g_p = f_p

glmn = so3.forward(g_before, g_p)

g = so3.inverse(glmn, g_p)

h, h_p = so3.convolve(f, f_p, g, g_p)
