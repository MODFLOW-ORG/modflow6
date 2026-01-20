import numpy as np

# power for Manning's hydraulic radius term
mpow = 2.0 / 3.0


def calculate_rectchan_mannings_discharge(
    conversion_factor, roughness, slope, width, depth
):
    """
    Calculate Manning's discharge for a rectangular channel.

    """
    area = width * depth
    return conversion_factor * area * depth**mpow * slope**0.5 / roughness


# n-point cross-section functions
def get_wetted_station(
    x0,
    x1,
    h0,
    h1,
    depth,
):
    """Get the wetted length in the x-direction"""
    # -- calculate the minimum and maximum depth
    hmin = min(h0, h1)
    hmax = max(h0, h1)

    # -- if depth is less than or equal to the minimum value the
    #    station length (xlen) is zero
    if depth <= hmin:
        x1 = x0
    # -- if depth is between hmin and hmax, station length is less
    #    than h1 - h0
    elif depth < hmax:
        xlen = x1 - x0
        dlen = h1 - h0
        if abs(dlen) > 0.0:
            slope = xlen / dlen
        else:
            slope = 0.0
        if h0 > h1:
            dx = (depth - h1) * slope
            xt = x1 + dx
            xt0 = xt
            xt1 = x1
        else:
            dx = (depth - h0) * slope
            xt = x0 + dx
            xt0 = x0
            xt1 = xt
        x0 = xt0
        x1 = xt1
    return x0, x1


def get_wetted_perimeter(
    x0,
    x1,
    h0,
    h1,
    depth,
):
    # -- calculate the minimum and maximum depth
    hmin = min(h0, h1)
    hmax = max(h0, h1)

    # -- calculate the wetted perimeter for the segment
    xlen = x1 - x0
    dlen = 0.0
    if xlen > 0.0:
        if depth >= hmax:
            dlen = hmax - hmin
        else:
            dlen = depth - hmin
            xlen = x1 - x0
    else:
        if depth > hmin:
            dlen = min(depth, hmax) - hmin
        else:
            dlen = 0.0

    return np.sqrt(xlen**2.0 + dlen**2.0)


def get_wetted_area(x0, x1, h0, h1, depth):
    # -- calculate the minimum and maximum depth
    hmin = min(h0, h1)
    hmax = max(h0, h1)

    # -- calculate the wetted area for the segment
    xlen = x1 - x0
    area = 0.0
    if xlen > 0.0:
        # -- add the area above hmax
        if depth > hmax:
            area = xlen * (depth - hmax)
        # -- add the area below zmax
        if hmax != hmin:
            if depth >= hmax:
                area += 0.5 * (hmax - hmin) * xlen
            elif depth > hmin:
                x0, x1 = get_wetted_station(x0, x1, h0, h1, depth)
                area += 0.5 * (depth - hmin) * (x1 - x0)
    return area


def wetted_area(
    x,
    h,
    depth,
    verbose=False,
):
    area = 0.0
    if x.shape[0] == 1:
        area = x[0] * depth
    else:
        for idx in range(0, x.shape[0] - 1):
            x0, x1 = x[idx], x[idx + 1]
            h0, h1 = h[idx], h[idx + 1]

            # get wetted area
            a = get_wetted_area(x0, x1, h0, h1, depth)
            area += a

            # write to screen
            if verbose:
                print(f"{idx}->{idx + 1} ({x0},{x1}) - perimeter={x1 - x0} - area={a}")

    return area


def add_wetted_vert(x, h, depth, vert_neighbs, idx):
    left_wet_len = 0
    right_wet_len = 0

    # left side
    if vert_neighbs[0]:
        idxm1 = idx - 1
        if h[idxm1] > depth:
            left_wet_len = depth - h[idx]
        else:
            left_wet_len = h[idxm1] - h[idx]

    # right side
    if vert_neighbs[1]:
        idxp1 = idx + 1
        idxp2 = idxp1 + 1
        if h[idxp2] > depth:
            right_wet_len = depth - h[idxp1]
        else:
            right_wet_len = h[idxp2] - h[idxp1]

    vert_len = left_wet_len + right_wet_len
    return vert_len


def wetted_perimeter(
    x,
    h,
    depth,
    verbose=False,
):
    perimeter = 0.0
    if x.shape[0] == 1:
        perimeter = x[0]
    else:
        for idx in range(0, x.shape[0] - 1):
            x0, x1 = x[idx], x[idx + 1]
            h0, h1 = h[idx], h[idx + 1]

            # get station data
            x0, x1 = get_wetted_station(x0, x1, h0, h1, depth)

            # get wetted perimeter
            perimeter += get_wetted_perimeter(x0, x1, h0, h1, depth)

            # set neighbor status
            vert_neighbs = is_neighb_vert(x, h, idx)

            # add wetted vertical neighbors if necessary
            if np.any(vert_neighbs):
                perimeter += add_wetted_vert(x, h, depth, vert_neighbs, idx)

            # write to screen
            if verbose:
                print(f"{idx}->{idx + 1} ({x0},{x1}) - perimeter={x1 - x0}")

    return perimeter


def is_neighb_vert(x, h, idx):
    cnt = len(x)

    # Assess left neighbor first
    if idx > 0:
        if cnt > 2:  # only x-sections w/ 3 or more pts may host a vertical side
            idxm1 = idx - 1
            if x[idxm1] == x[idx] and h[idxm1] != h[idx]:
                leftvert = True
            else:
                leftvert = False
        else:
            leftvert = False
    else:
        leftvert = False

    # Assess right neighbor
    idxp1 = idx + 1
    idxp2 = idxp1 + 1
    if cnt > idxp2:
        if x[idxp1] == x[idxp2] and h[idxp1] != idxp2:
            rightvert = True
        else:
            rightvert = False
    else:
        rightvert = False

    return (leftvert, rightvert)


def manningsq(
    x,
    h,
    depth,
    roughness=0.01,
    slope=0.001,
    conv=1.0,
):
    if isinstance(roughness, float):
        roughness = np.ones(x.shape, dtype=float) * roughness
    if x.shape[0] > 1:
        q = 0.0
        for i0 in range(x.shape[0] - 1):
            i1 = i0 + 1

            # get station data
            x0, x1 = get_wetted_station(x[i0], x[i1], h[i0], h[i1], depth)

            perimeter = get_wetted_perimeter(x0, x1, h[i0], h[i1], depth)

            # set neighbor status
            vert_neighbs = is_neighb_vert(x, h, i0)

            # add wetted vertical neighbors if necessary
            if np.any(vert_neighbs):
                perimeter += add_wetted_vert(x, h, depth, vert_neighbs, i0)

            area = get_wetted_area(x[i0], x[i1], h[i0], h[i1], depth)
            if perimeter > 0.0:
                radius = area / perimeter
                q += conv * area * radius**mpow * slope**0.5 / roughness[i0]
    else:
        perimeter = wetted_perimeter(x, h, depth)
        area = wetted_area(x, h, depth)
        radius = 0.0
        if perimeter > 0.0:
            radius = area / perimeter
        q = conv * area * radius**mpow * slope**0.5 / roughness[0]
    return q


def get_depths(
    flows,
    x,
    h,
    roughness=0.01,
    slope=0.001,
    conv=1.0,
    dd=1e-4,
    verbose=False,
):
    if isinstance(flows, float):
        flows = np.array([flows], dtype=float)
    if isinstance(roughness, float):
        roughness = np.ones(x.shape, dtype=float) * roughness
    depths = np.zeros(flows.shape, dtype=float)
    for idx, q in enumerate(flows):
        depths[idx] = qtodepth(
            x, h, q, roughness=roughness, slope=slope, conv=conv, dd=dd, verbose=False
        )

    return depths


def qtodepth(
    x,
    h,
    q,
    roughness=0.01,
    slope=0.001,
    conv=1.0,
    dd=1e-4,
    verbose=False,
):
    h0 = 0.0
    q0 = manningsq(x, h, h0, roughness=roughness, slope=slope, conv=conv)
    r = q0 - q

    iter = 0
    if verbose:
        print(f"iteration {iter:>2d} - residual={r}")
    while abs(r) > 1e-12:
        q1 = manningsq(x, h, h0 + dd, roughness=roughness, slope=slope, conv=conv)
        dq = q1 - q0
        if dq != 0.0:
            derv = dd / (q1 - q0)
        else:
            derv = 0.0
        h0 -= derv * r
        q0 = manningsq(x, h, h0, roughness=roughness, slope=slope, conv=conv)
        r = q0 - q

        iter += 1
        if verbose:
            print(f"iteration {iter:>2d} - residual={r}")
        if iter > 100:
            break
    return h0
