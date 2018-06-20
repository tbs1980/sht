# sht
Portable header only shperical harmonic transforms.

## build
```commandline
mkdir build; cd build
cmake -GNinja ..\
ninja
```

## `alm2map` or synthesis

```
for m = [0:mmax]
  for l = [m:lmax]
    compute_alm_and_blm() # Eqn 14
  for r = [0:rmax]
    for l = [m:lmax]
      compute_s_lambda_lm(theta_r) # Eqn 13
      accumulate_F(m,theta_y) # Eqn 8

for r = [0:rmax]
  compute_map_using_fft(x,y)
```
