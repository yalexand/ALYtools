using Images

################################################
function radon(A::AbstractMatrix, angles)
    angles *= 1/360*2*pi
    w, h = size(A)
    d = floor(Int, sqrt(w*w+h*h)/2)  # half-length of the diagonal
    Nr = 2d+1
    R = linspace(-Nr/2, Nr/2, Nr)

    sinogram = zeros(Nr, length(angles))

    for j = 1:length(angles)
        a = angles[j]
        θ = a + pi/2
        sinθ, cosθ = sin(θ), cos(θ)
        for k = 1:length(R)
            rk = R[k]
            x0, y0 = rk*cosθ + w/2, rk*sinθ + h/2
            # for a line perpendicular to this displacement from the center,
            # determine the length in either direction within the boundary of the image
            lfirst, llast = interior(x0, y0, sinθ, cosθ, w, h, R)
            tmp = 0.0
            @inbounds for l in lfirst:llast
                rl = R[l]
                x, y = x0 - rl*sinθ, y0 + rl*cosθ
                ix, iy = trunc(Int, x), trunc(Int, y)
                fx, fy = x-ix, y-iy
                tmp += (1-fx)*((1-fy)*A[ix,iy] + fy*A[ix,iy+1]) +
                       fx*((1-fy)*A[ix+1,iy] + fy*A[ix+1,iy+1])
            end
            sinogram[k,j] = tmp
        end
    end
    return sinogram
end

########################################
function interior(x0, y0, sinθ, cosθ, w, h, R::Range)
    rx1, rxw = (x0-1)/sinθ, (x0-w)/sinθ
    ry1, ryh = (1-y0)/cosθ, (h-y0)/cosθ
    rxlo, rxhi = minmax(rx1, rxw)
    rylo, ryhi = minmax(ry1, ryh)
    rfirst = max(minimum(R), ceil(Int,  max(rxlo, rylo)))
    rlast  = min(maximum(R), floor(Int, min(rxhi, ryhi)))
    Rstart, Rstep = first(R), step(R)
    lfirst, llast = ceil(Int, (rfirst-Rstart)/Rstep), floor(Int, (rlast-Rstart)/Rstep)
    if lfirst == 0
        lfirst = length(R)+1
    end
    # Because of roundoff error, we still can't trust that the forward
    # calculation will be correct. Make adjustments as needed.
    while lfirst <= llast
        rl = R[lfirst]
        x, y = x0 - rl*sinθ, y0 + rl*cosθ
        if (1 <= x < w) && (1 <= y < h)
            break
        end
        lfirst += 1
    end
    while lfirst <= llast
        rl = R[llast]
        x, y = x0 - rl*sinθ, y0 + rl*cosθ
        if (1 <= x < w) && (1 <= y < h)
            break
        end
        llast -= 1
    end
    lfirst, llast
end

#########################################
function iradon(A,angles) # input angles [degree]

  angles *= 1/360*2*pi

  angles = angles + pi/2

  N, nAngles = size(A)

  I = zeros(N,N); # reconstruction

  x = linspace(-0.5,0.5,N)

  filter = abs(linspace(-1, 1, N))

  # FT domain filtering
  for t=1:length(angles)
      #fhat = fftshift(fft(slice(A,:,t)))
      fhat = fftshift(fft(A[:,t]))
      A[:,t] = real(ifft(ifftshift(fhat.*filter)))
  end

  XCOS = zeros(N,length(angles))
  XSIN = zeros(N,length(angles))
  for k=1:N
    for a=1:length(angles)
      XCOS[k,a]=x[k]*cos(angles[a])
      XSIN[k,a]=x[k]*sin(angles[a])
    end
  end

  @inbounds for t=1:nAngles
          for n=1:N
              xs = XSIN[n,t]
              for m=1:N
                  r = XCOS[m,t] + xs
                  index = Int64(round((r/sqrt(2.0) + 0.5)*N)) # sqrt(2) magnified
                  if index>0 && index<=N
                      I[m,n] += A[index,t]
                  end
              end
          end
      end

  N_ = Int64(round(N/sqrt(2.)))
  I = imresize(I,N_,N_)

  return I
end
