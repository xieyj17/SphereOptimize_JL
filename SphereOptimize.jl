@everywhere function regularized_Theta(ang)
    ang = ang % (2*pi)
    while ang < 0
        ang += 2*pi
    end
    ang = ang % (2*pi)
    return ang
end

@everywhere function med_from_Sphere(theta, i)
    res = 1
    if i <= length(theta)
        temp_theta = theta[1:i]
        if i == 1
            res = cos(temp_theta[1])
        else
            for p in 1:(i-1)
                res = res * sin(temp_theta[p])
            end
            res = res * cos(temp_theta[i])
        end
    else
        temp_theta = theta
        for p in 1:(i-2)
            res *= sin(temp_theta[p])
        end
        res = res * cos(temp_theta[i-1])
    end
    return res
end


@everywhere function from_Sphere(theta)
    s = zeros(length(theta)+1)
    for k in 1:length(theta)
        theta[k] = regularized_Theta(theta[k])
    end
    for i in 1:length(s)
        s[i] = med_from_Sphere(theta, i)
    end
    return s
end


@everywhere function to_Sphere(s)
    theta = zeros(length(s)-1)
    theta[1] = acos(s[1])
    for i in 2:length(theta)
        theta[i] = acos(s[i] / sqrt(sum([w^2 for w in s[i:length(s)]])))
    end
    if s[length(s)] < 0
        theta[length(theta)] = 2*pi - theta[length(theta)]
    end
    if ismissing(theta)
        theta[[ismissing(x) for x in theta]] = 0
    end
    return theta
end

function Sphere_Optimize(par, fn, neighbor = NaN, ...)
    function temp_fn(t)
        s = from_Sphere(t)
        res = fn(s)
        return res
    end

    theta = to_Sphere(par)
    k = optimize(temp_fn, theta, LBFGS())
    ntheta = Optim.minimizer(k)
    ntheta = [regularized_Theta(t) for t in ntheta]
    par = to_Sphere(ntheta)
    nv = Optim.minimum(k)

    return [nv, par]
end
