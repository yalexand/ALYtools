function hR = hR2(x, angles, angles2)

hR = [radon(x, angles)  radon(fliplr(flipud(x)),angles2)];