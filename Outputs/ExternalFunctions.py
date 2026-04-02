def e_geo(traces, x, y):
    """
    Calculate the geomagnetic component from the electric field in the shower plane,
    i.e. the electric field should be in the (vxB, vxvxB, v) CS

    Parameters
    ----------
    traces : np.ndarray
        The traces in the shower plane, shaped as (samples, polarisations)
    x : float
        The antenna position along the vxB axis
    y : float
        The antenna position along the vxvxB axis

    Returns
    -------
    e_geo : np.ndarray
        The geomagnetic component of the electric field
    """
    return traces[:, 1] * x / y - traces[:, 0]


def e_ce(traces, x, y):
    """
    Calculate the charge-excess (or Askaryan) component of electric field in the shower plane,
    i.e. the electric field should be in the (vxB, vxvxB, v) CS

    Parameters
    ----------
    traces : np.ndarray
        The traces in the shower plane, shaped as (samples, polarisations)
    x : float
        The antenna position along the vxB axis
    y : float
        The antenna position along the vxvxB axis

    Returns
    -------
    e_ce : np.ndarray
        The charge-excess component of the electric field
    """
    return -traces[:, 1] * np.sqrt(x**2 + y**2) / y


def e_to_geo_ce(traces, x, y):
    """
    Decouples the electric field in the shower plane, i.e. the electric field should be in the (vxB, vxvxB, v) CS,
    into the geomagnetic and charge-excess components.

    Parameters
    ----------
    traces : np.ndarray
        The traces in the shower plane, shaped as (samples, polarisations)
    x : float
        The antenna position along the vxB axis
    y : float
        The antenna position along the vxvxB axis

    Returns
    -------
    e_geo : np.ndarray
        The geomagnetic component of the electric field
    e_ce : np.ndarray
        The charge-excess component of the electric field
    """
    return e_geo(traces, x, y), e_ce(traces, x, y)


def geo_ce_to_e(my_e_geo, my_e_ce, x, y):
    """
    Convert the geomagnetic and charge-excess components to a three-dimensional electric field in
    the shower plane, i.e. the (vxB, vxvxB, v) CS.

    Note that the v-component is simply set to zero.

    Parameters
    ----------
    my_e_geo: np.ndarray
        The electric field traces of the geomagnetic component, shaped as (n_antennas, n_samples, n_slices,)
    my_e_ce: np.ndarray
        The electric field traces of the charge-excess component, shaped as (n_antennas, n_samples, n_slices,)
    x : float
        The antenna position along the vxB axis, shaped as (n_antennas, )
    y : float
        The antenna position along the vxvxB axis, shaped as (n_antennas, )

    Returns
    -------
    e_field : np.ndarray
        The three-dimensional electric field, in the shower plane CS, shaped as (samples, polarisations)
    """
    if hasattr(x, "__len__"):
        x = np.asarray(x)
        y = np.asarray(y)
    else:
        x = np.asarray([x])
        y = np.asarray([y])

    # expand the x and y positions to match the traces
    expand_axes = (1, 2) if my_e_geo.ndim == 3 else (1,)
    x_arr = np.expand_dims(x, axis=expand_axes)
    y_arr = np.expand_dims(y, axis=expand_axes)

    trace_vB = -1 * (
        my_e_geo
        + my_e_ce
        * x_arr
        / np.sqrt(x_arr ** 2 + y_arr ** 2)
    )
    trace_vvB = (
        -1
        * my_e_ce
        * y_arr
        / np.sqrt(x_arr ** 2 + y_arr ** 2)
    )
    trace_v = np.zeros_like(trace_vB)

    return np.stack((trace_vB, trace_vvB, trace_v), axis=-1)