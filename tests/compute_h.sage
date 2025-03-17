from sage.rings.power_series_ring_element import PowerSeries

n = 40

R.<z> = PowerSeriesRing(ZZ, 'z')
p = PowerSeries.polynomial((1 + z^2)^(n+1)/(1-z)^n)