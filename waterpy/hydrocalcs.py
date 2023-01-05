"""
Module of hydrologic calculations.

References:

Engineering and Design - Runoff from Snowmelt
U.S. Army Corps of Engineers
Engineering Manual 1110-2-1406
https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/snow/COEemSnowmeltRunoff.pdf

Estimating Runoff
U.S. Department of Agriculture - Urban Hydrology for Small Watersheds
Technical Release 55 (TR55), June 1986
https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1044171.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats


def pet(dates, temperatures, latitude, calib_coeff, method="hamon"):
    """Calculate potential evapotranspiration for various methods.

    :param dates: An array of python datetimes
    :type dates: numpy.ndarray
    :param temperatures: An array of temperatures , in degrees Celsius
    :type temperatures: numpy.ndarray
    :param latitude: A latitude, in decimal degrees
    :type latitude: float
    :param calib_coeff: Calibration coefficient (KPEC), dimensionless
    :type calib_coeff: float
    :return pet: array of pet values, in millimeters per day
    :rtype pet: numpy.ndarray
    """
    if method.lower() == "hamon":
        return pet_hamon(dates, temperatures, latitude, calib_coeff)


def pet_hamon(dates, temperatures, latitude, calib_coeff=1.2):
    """Calculate the amount of potential evapotranspiration in millimeters
    per day using the Hamon equation.

    :param dates: An array of python datetimes
    :type dates: numpy.ndarray
    :param temps: An array of temps, in degrees Celsius
    :type temps: numpy.ndarray
    :param latitude: A latitude, in decimal degrees
    :type latitude: float
    :param calib_coeff: Calibration coefficient (KPEC), dimensionless
    :type calib_coeff: float
    :return pet: array of pet values, in millimeters per day
    :rtype pet: numpy.ndarray

    .. notes::
       Equation:

       (1) PET = 0.1651 * Ld * RHOSAT * KPEC

       where

           PET          daily potential evapotranspiration (mm/day)
           Ld           daytime length (daylight hours), time from sunrise
                        to sunset in multiples of 12 hours (hours)
           RHOSAT       saturated vapor density at the daily mean air
                        temperature (T), g/m**3
           KPEC         calibration coefficient, dimensionless
                        set to 1.2 for Southeastern United States

       Sub-equations:

       (2) RHOSAT = 216.7 * ESAT / (T + 273.3)

       (3) ESAT = 6.108 * EXP(17.26939 * T / (T + 237.3))

       where

           T       daily mean air temperature, (celsius)
           ESAT    saturated vapor pressure at the given T, (mb)

       (4) Ld = (w / 15) * 2

       where

           w       the sunset hour angle (degrees); Earth moves 15° per hour

       (5) w = arccos(-1 * tan(declination) * tan(latitude))

       where

           latitude           angle distance of a place north  or south of the
                              earth's equator (radians)

           declination        angle between the Sun's rays and the equatorial
                              plane (radians); the declination of the Earth is
                              the angular distance at solar noon between the
                              Sun and the Equator, north-positive

       (6) declination = 23.45° * sin((360 * (284 + N) / 365))

       where

           N                  number of days after January 1 (the Julian Day)

       References:
           - Lu et al. (2005). A comparison of six potential evapotranspiration
             methods for reginal use in the southeastern United States.

           - Brock, Thomas. (1981). Calculating solar radiation for ecological
             studies. Ecological Modeling.

           - https://en.wikipedia.org/wiki/Position_of_the_Sun
    """
    if len(dates) != len(temperatures):
        raise IndexError(
            "Length of dates: {}\n",
            "Length of temperatures: {}\n",
            "Lengths of dates and temperatures must be equal"
            "".format(len(dates), len(temperatures))
        )

    DEG2RAD = np.pi/180
    RAD2DEG = 180/np.pi

    pet = []
    for date, temperature in zip(dates, temperatures):

        # Declination
        # Note: using python datetimes whih have .timetuple() method
        day_num = date.timetuple().tm_yday
        angle = 360 * ((284 + day_num) / 365) * DEG2RAD
        declination = (23.45 * DEG2RAD) * np.sin(angle)

        # calculate sunset hour angle in degrees (w)
        sunset_hour_angle = (
            np.arccos(-1 * np.tan(declination) * np.tan(latitude * DEG2RAD))
            * RAD2DEG
        )

        # calculate daytime length in 12 hour unit (Ld)
        daytime_length = abs((sunset_hour_angle / 15) * 2) / 12

        # calculate saturated vapor pressure (ESAT)
        saturated_vapor_pressure = (
            6.108 * np.exp((17.26939 * temperature) / (temperature + 237.3))
        )

        # calculate saturated vapor density (RHOSAT)
        saturated_vapor_density = (
            (216.7 * saturated_vapor_pressure) / (temperature + 273.3)
        )

        if temperature <= 0:
        # If temperatures are below freezing PET shouldn't exist.

            potential_evapotranspiration = 0

        else:
        # calculate potential evapotranspiration
            potential_evapotranspiration = (
                0.1651 * daytime_length * saturated_vapor_density * calib_coeff
            )

        # add to list
        pet.append(potential_evapotranspiration)

    # convert list to numpy array
    pet = np.array(pet)

    return pet


def snowmelt(precipitation,
             temperatures,
             temperature_cutoff,
             snowmelt_rate_coeff_with_rain,
             snowmelt_rate_coeff,
             timestep_daily_fraction):
    """Snow melt routine.

    Uses snow water equivalence which is the amount of water contained
    within the snowpack. Snow water equivalence is assumed to be 10% of
    the snow water density.

    In addition to melting from the "top" of the snow pack due to temperature,
    it is assumed that some melting occurs at the "bottom" of the snow pack
    each day regardless of temperature due to geothermal energy as mentioned
    from:

    Heat conducted from the ground becomes an energy source, which in
    turn causes snow melt from bottom of snowpack. Male and Gray estimated melt
    rates on the order of 0.025 to 0.076 cm/day. (0.25 -0.76 mm)

    Dingman estimated the ground melt parameter default ~ 0.35 mm

    Dingman S.L.(2002) Physical Hydrology Printace Hall 2nd ed.

    Male D.H. and Gray D.M.(1981) "Snowcover Ablation and Runoff" in Handbook
    of Snow" Principles Processes, Management and Use, D.M. Gray and D.H. Male,
    ed., Pergamon Press:360-436

    :param precipitation: Precipitation rates, in millimeters per day
    :type precipitation: numpy.ndarray
    :param temperatures: Temperatures, in degrees Fahrenheit
    :type temperatures: numpy.ndarray
    :param temperature_cutoff: Temperature when melt begins,
                               in degrees Fahrenheit
    :type temperature_cutoff: float
    :param snowmelt_rate_coeff_with_rain: Snowmelt coefficient when raining,
                                          1/degrees Fahrenheit
    :type snowmelt_rate_coeff_with_rain: float
    :param snowmelt_rate_coeff: Snowmelt rate coefficient (often variable),
                                in inches per degree Fahrenheit
    :type snowmelt_rate_coeff: float
    :param timestep_daily_fraction: Model timestep as a fraction of a day
    :type timestep_daily_fraction: float
    :return: Tuple of arrays of adjusted precipitation, snowmelt,
             snowpack values, and snow water equivalence, each array
             is in millimeters per day
    :rtype: Tuple

    """
    # Snow water equivalence is assumed to be 10% of the snow water density.
    snow_water_equivalence_factor = 0.1

    precip_inches = precipitation / 25.4  # mm to inches

    snowprecip = []
    snowmelts = []
    snowpacks = []

    snowmelt = 0
    snowpack = 0
    for temp, precip_inch in zip(temperatures, precip_inches):

        # If temp is high enough then there is snowmelt,
        # calculate amount of snowmelt, else no snowmelt,
        # snowpack accumulates by full precip amount
        if temp >= temperature_cutoff:
            # If it is raining, calculate snowmelt with rain,
            # else calculate snowmelt without rain
            if precip_inch > 0:
                snowmelt = snowmelt_rain_on_snow_heavily_forested(
                    precip_inch,
                    temp,
                    temperature_cutoff,
                    snowmelt_rate_coeff_with_rain
                )

                # adjust daily snowmelt to same timestep as precip and temp
                snowmelt = snowmelt * timestep_daily_fraction

            else:
                snowmelt = snowmelt_temperature_index(
                    temp,
                    temperature_cutoff,
                    snowmelt_rate_coeff
                )

                # adjust daily snowmelt to same timestep as precip and temp
                snowmelt = snowmelt * timestep_daily_fraction

            # If there is more snowmelt than snowpack available to melt, then
            # add the amount of snowpack available to the snowmelt
            if snowmelt >= snowpack:
                snowmelt = snowpack

            # Remove the amount of snowmelt from the snowpack,
            # and add the amount of snowmelt to the current precip amount
            snowpack = snowpack - snowmelt
            precip_inch = precip_inch + (
                snowmelt * snow_water_equivalence_factor
            )

        # If temp is too cold for melting, then add the precip (snow) amount
        # to the snow pack and assume no water infiltrates on cold days by
        # setting precip to zero
        else:
            snowpack = snowpack + (precip_inch / snow_water_equivalence_factor)
            precip_inch = 0

        # Apply geothermal melting from "bottom" of snow pack
        # Note: 0.51 mm/day (or 0.02 inches/day) is the midway point between
        # range of estimated melt rates from Male and Gray (see reference
        # in docstring). Male and Gray estimated melt rates on the order of
        # 0.025 to 0.076 cm/day or
        # 0.25 to 0.76 mm/day or
        # 0.01 to 0.03 inches/day
        melt_rate_male_gray = 0.02  # inches/day
        if snowpack >= melt_rate_male_gray:
            snowpack = snowpack - melt_rate_male_gray
            precip_inch = precip_inch + (
                melt_rate_male_gray * snow_water_equivalence_factor
            )
        else:
            precip_inch = precip_inch + (
                snowpack * snow_water_equivalence_factor
            )
            snowpack = 0

        snowprecip.append(precip_inch)
        snowmelts.append(snowmelt)
        snowpacks.append(snowpack)

    snowprecip = np.array(snowprecip) * 25.4  # inches to mm
    snowmelts = np.array(snowmelts) * 25.4  # inches to mm
    snowpacks = np.array(snowpacks) * 25.4  # inches to mm

    snow_water_equivalence = snowpacks * snow_water_equivalence_factor

    return snowprecip, snowmelts, snowpacks, snow_water_equivalence


def snowmelt_rain_on_snow_heavily_forested(precipitation,
                                           temperatures,
                                           temperature_cutoff=32.0,
                                           rain_melt_coeff=0.007):
    """Calculate the amount of snowmelt rain-on-snow situations in
    heavily forested areas using a generalized equation for rain-on-snow
    situations in heavily forested areas (the mean canopy cover is greater
    than 80%). This snowmelt calculation is from the family of energy budget
    solutions.

    :param precipitation: Precipitation rates, in inches/day
    :type precipitation: numpy.ndarray
    :param temperatures: Temperatures, in degrees Fahrenheit
    :type temperatures: numpy.ndarray
    :param temperature_cutoff: Temperature when melt begins,
                               in degrees Fahrenheit
    :type temperature_cutoff: float
    :param rain_melt_coeff: Snowmelt coefficient when raining,
                             1/degrees Fahrenheit
    :type rain_melt_coeff: float
    :return snowmelt: Snowmelt values, in inches per day
    :rtype snowmelt: numpy.ndarray

    .. note::

    Equation:

    M = (0.074 + 0.007 * P_r) * (T_a - 32) + 0.05

    where

    M           snowmelt, inches/day
    P_r         rate of precipitation, inches/day
    T_a         temperature of saturated air, at 3 meters (10 ft) level,
                degrees Fahrenheit

    Reference:

    Engineering and Design - Runoff from Snowmelt
    U.S. Army Corps of Engineers
    Engineering Manual 1110-2-1406
    Chapter 5-3. Generalized Equations, Rain-on-Snow Situations, Equation 5-20
    https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/snow/COEemSnowmeltRunoff.pdf
    """
    return (
        (0.074 + rain_melt_coeff * precipitation)
        * (temperatures - temperature_cutoff) + 0.05
    )


def snowmelt_rain_on_snow_open_to_partly_forested(precipitation,
                                                  temperatures,
                                                  winds,
                                                  temperature_cutoff=32.0,
                                                  rain_melt_coeff=0.007,
                                                  basin_wind_coeff=0.5):
    """Calculate the amount of snowmelt rain-on-snow situations in
    partially forested areas using a generalized equation for
    rain-on-snow situations in partially forested areas (the mean canopy
    cover is greater than 10%-80%). Snowmelt calculation is from the family
    of energy budget solutions.

    :param precipitation: Precipitation rates, in inches/day
    :type precipitation: numpy.ndarray
    :param temperatures: Temperatures of saturated air,
                         in degrees Fahrenheit
    :type temperatures: numpy.ndarray
    :param winds: Winds, in miles per hour
    :type winds: numpy.ndarray
    :param temperature_cutoff: Temperature when melt begins,
                               in degrees Fahrenheit
    :type temperature_cutoff: float
    :param rain_melt_coeff: Snowmelt coefficient when raining,
                             1/degrees Fahrenheit
    :type rain_melt_coeff: float
    :param basin_wind_coeff: Basin wind exposure coefficient, fraction
    :type basin_wind_coeff: float
    :return snowmelt: Snowmelt values, in inches per day
    :rtype snowmelt: numpy.ndarray

    .. note::

    Equation:

    M = (0.029 + 0.0084 * k * v + 0.007 * P_r) * (T_a - 32) + 0.09

    where

    M           snowmelt, inches/day
    k           basin wind exposure coefficient, unitless
    v           wind velocity, miles/hour
    P_r         rate of precipitation, inches/day
    T_a         temperature of saturated air, at 3 meters (10 ft) level,
                degrees Fahrenheit

    Reference:

    Engineering and Design - Runoff from Snowmelt
    U.S. Army Corps of Engineers
    Engineering Manual 1110-2-1406
    Chapter 5-3. Generalized Equations, Rain-on-Snow Situations, Equation 5-19
    https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/snow/COEemSnowmeltRunoff.pdf
    """
    return (
        (0.029
         + (0.0084 * basin_wind_coeff * winds)
         + (rain_melt_coeff * precipitation))
        * (temperatures - temperature_cutoff) + 0.09
    )


def snowmelt_temperature_index(temperatures,
                               temperature_cutoff=32.0,
                               melt_rate_coeff=0.06):
    """Calculate the amount of snowmelt using a temperature index method,
    also called degree-day method. This method has its limitations as noted
    in the reference.

    :param temperatures: Temperatures, in degrees Fahrenheit
    :type temperatures: numpy.ndarray
    :param temperature_cutoff: Temperature when melt begins,
                               in degrees Fahrenheit
    :type temperature_cutoff: float
    :param melt_rate_coeff: Snowmelt rate coefficient (often variable),
                            in inches per degree Fahrenheit
    :type melt_rate_coeff: float
    :return snowmelt: Snowmelt values, in inches per day
    :rtype snowmelt: numpy.ndarray

    .. note::

    Equation:

    M = C_m * (T_a - T_b)

    where

    M           snowmelt, inches per period
    C_m         melt-rate coefficient, inches/(degree Fahrenheit/period)
    T_a         air temperature, degrees Fahrenheit
    T_b         base temperature, degrees Fahrenheit

    Reference:

    Engineering and Design - Runoff from Snowmelt
    U.S. Army Corps of Engineers
    Engineering Manual 1110-2-1406
    Chapter 6-1. Generalized Equations, Rain-on-Snow Situations, Equation 6-1
    https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/snow/COEemSnowmeltRunoff.pdf
    """
    return melt_rate_coeff * (temperatures - temperature_cutoff)


def runoff(grow_season, precipitation, curve_number, amc):
    """Calculate the amount of runoff using the SCS runoff curve number method.
    :param temperature: temperature, in C
    :type temperature: numpy.ndarray
    :param precipitation: precipitation, in mm
    :type precipitation: numpy.ndarray
    :param curve_number: curve number, dimensionless
    :type precipitation: float
    :return runoff: runoff, in mm
    :rtype runoff: numpy.ndarray

    .. note::

    Equation:

    Q = (P - I_a)**2 / ((P - I_a) + S)

    where

    Q           runoff, inches
    P           rainfall, inches
    S           potential maximum retention after runoff begins, inches
    I_a         initial abstraction, inches

    using approximation

    I_a = 0.2 * S

    and

    S = (1000 / CN) - 10

    where

    CN          curve number [0-100], dimensionless

    Refactored Equation:

    Q = (P - 0.2*S)**2 / (P + 0.8*S)


    Reference:

    Estimating Runoff
    U.S. Department of Agriculture - Urban Hydrology for Small Watersheds
    Technical Release 55 (TR55), June 1986
    https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1044171.pdf
    """
    # might need a 0.1 delay timing on the storage in this section.

    if grow_season is True:
        # non-growing season
        if amc > 12.7:
            if amc > 27.94:
                # wet and cold
                curve_number = 23 * curve_number / (10 + 0.13 * curve_number)
            else:
                # normal and cold
                curve_number = curve_number
        else:
            # dry and cold
            curve_number = curve_number * 4.2 / (10 - 0.058 * curve_number)
    else:
        if amc > 35.56:
            if amc > 53.34:
                # wet and warm
                curve_number = 23 * curve_number / (10 + 0.13 * curve_number)
            else:
                # normal and warm
                curve_number = curve_number
        else:
            # dry and warm
            curve_number = curve_number * 4.2 / (10 - 0.058 * curve_number)

    precip_inches = precipitation / 25.4  # mm to inches

    potential_retention = (1000 / curve_number) - 10
    runoff_inches = (
        (precip_inches - 0.2 * potential_retention)**2
        / (precip_inches + 0.8 * potential_retention)
    )
    runoff = runoff_inches * 25.4  # inches to mm

    return runoff

def weighted_mean(values, weights):
    """Calculate the weighted mean.

    :param values: Array of values
    :type values: numpy.ndarray
    :param weights: Array of weights
    :type weights: numpy.ndarray
    :rtype: float
    """
    weighted_mean = (values * weights).sum() / weights.sum()

    return weighted_mean


def absolute_error(observed, modeled):
    """Calculate the absolute error between two arrays.

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: numpy.ndarray
    """
    error = observed - modeled

    return error


def mean_squared_error(observed, modeled):
    """Calculate the mean square error between two arrays.

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: float
    """
    error = absolute_error(observed, modeled)
    mse = np.mean(error**2)

    return mse


def relative_error(observed, modeled):
    """Calculate the relative change between two arrays.

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: numpy.ndarray
    """
    error = absolute_error(observed, modeled) / observed

    return error


def percent_error(observed, modeled):
    """Calculate the percent error between two arrays.

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: numpy.ndarray
    """
    error = relative_error(observed, modeled) * 100

    return error


def percent_difference(observed, modeled):
    """Calculate the percent difference between two arrays.

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: numpy.ndarray
    """
    mean = np.mean((observed, modeled), axis=0)
    percent_diff = ((observed - modeled) / mean) * 100

    return percent_diff


def r_squared(observed, modeled):
    """Calculate the Coefficient of Determination. Used to indicate how well
    data points fit a line or curve. Use numpy.coeff for computation.

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: float
    """
    r = np.corrcoef(observed, modeled)[0, 1]
    coefficient = r**2

    return coefficient


def nash_sutcliffe(observed, modeled):
    """Calculate the Nash-Sutcliffe (model efficiency coefficient).
    Used to assess the predictive power of hydrological models.

    E = 1 - sum((observed - modeled) ** 2)) / (sum((observed - mean_observed)**2 )))

    :param observed: Array of observed data
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :type modeled: numpy.ndarray
    :rtype: float
    """
    mean_observed = np.mean(observed)
    numerator = np.sum((observed - modeled) ** 2)
    denominator = np.sum((observed - mean_observed)**2)
    coefficient = 1 - (numerator/denominator)

    return coefficient


def flow_duration(values):
    """Calculate the exceedance probabilities for a set of values for use in
    plotting a flow duration curve.

    :param observed: Array of flow values
    :type observed: numpy.ndarray
    :param modeled: Array of modeled data
    :return tuple: Tuple of probabilities, sorted values
    :rtype: tuple
    """
    # Sort the values
    values_sorted = np.sort(values)

    # Rank data from smallest to largest
    ranks = stats.rankdata(values_sorted, method="average")

    # Reverse the order
    ranks = ranks[::-1]

    # Compute the exceedance probabilities
    probabilities = (
        [(ranks[i] / (len(values) + 1)) * 100 for i in range(len(values))]
    )

    return probabilities, values_sorted


def randomize(value, size=24):
    """Generate pseudo random values of a certain size that sum to 1.0.
    Used to create 24 hourly random values from a single daily value.

    :param values: Value
    :type values: float
    :param size: Size of array wanted to be returned, default is 24 (hours)
    :type size: int
    :return randomized: Array of random values.
    :rtype: numpy.ndarray

    """
    # Set the seed to reproduce results
    np.random.seed(1)

    # Create random array of a certain size
    random_array = np.random.random(size)

    # Normalize to make the sum of the array equal to 1.0
    random_array_normalized = random_array / random_array.sum()

    # Distribute value across array
    randomized = random_array_normalized * value

    return randomized

def create_rain_array(rain_file, values):
    """Creates a random array of rainfall values for hourly distribution
    :param rain_file: csv file
    :type rain_file: hourly rainfall distribution file defined in config file.
    :param values: Value
    :type values: np.ndarray
    :return: rain_dist
    :rtype: np.ndarray
    """

    np.random.seed(10)
    df = pd.read_csv(rain_file)
    val_len = len(values) + 1
    val_array = np.array(range(0, val_len))

    rain_dist = []
    high = len(df.columns) - 1

    rand_array = np.random.randint(low=0, high=high, size=len(val_array))

    for day in rand_array:
        day_dist = df[df.columns[day]].values
        for value in day_dist:
            rain_dist.append(value)

    return np.array(rain_dist)








def randomize_daily_to_hourly(values):
    """Generate pseudo random values of a certain size that sum to 1.0.
    Used to create hourly random values from daily values.

    :param values: Daily values to randomize into hourly values.
    :type values: numpy.ndarray
    :rtype: numpy.ndarray
    """
    # Set the seed to reproduce results
    np.random.seed(1)

    # Create an array of random hourly values from daily values
    randomized = []
    for value in values:
        # Create random array of 24 hours and add some zeros
        n_zeros = 0 #np.random.randint(1, 12)
        n_non_zeros = 24 - n_zeros

        rand = np.random.random(n_non_zeros)
        rand = np.append(rand, np.zeros(n_zeros))

        # Shuffle the array to that the zeros are dispersed throughout the day
        np.random.shuffle(rand)

        # Normalize to make the sum of the array equal to 1.0
        rand_normalized = rand / rand.sum()

        # Distribute value across array
        rand_with_value = rand_normalized * value

        # Loop through random array with value and append random values to list
        for rand_value in rand_with_value:
            randomized.append(rand_value)

    return np.array(randomized)


def copy_daily_to_hourly(values):
    """Create hourly values from daily values by copying each daily value
    and spreading it across and hourly array.

    :param values: Daily values to copy into hourly values.
    :type values: numpy.ndarray
    :rtype: numpy.ndarray
    """
    # Create an array of hourly values from daily values
    copied = []
    for value in values:
        # Create array for 24 hours and distribute (copy) daily value to hourly
        hourly_array = np.ones(24) * value

        # Loop through array and append hourly values to list
        for hourly_value in hourly_array:
            copied.append(hourly_value)

    return np.array(copied)


def chop_daily_to_hourly(values):
    """Create hourly values from daily values by evenly chopping each daily value
    and spreading it across an hourly array.  This is akin to the fortran version
    of the model.

    :param values: Daily values from timeseries file.
    :type values: numpy.ndarray
    :rtype: numpy.ndarray
    """
    copied = []
    for value in values:
        hourly_array = np.ones(24) * value / 24
        for hourly_value in hourly_array:
            copied.append(hourly_value)

    return np.array(copied)

def chop_daily_to_hourly_precip(values, rain_array):
    """Create hourly values from daily values by evenly chopping each daily value
    and spreading it across an hourly array.  This is akin to the fortran version
    of the model.

    :param values: Daily values from timeseries file.
    :type values: numpy.ndarray
    :rtype: numpy.ndarray
    """
    copied = []
    i = 0
    for value in values:
        hourly_array = np.ones(24) * value * rain_array[i]
        i += 1
        for hourly_value in hourly_array:
            copied.append(hourly_value)

    return np.array(copied)

def sum_hourly_to_daily(values, minmax=False):
    """Create daily values from hourly values by summing in every 24th element
    in the array.

    :param values: Hourly values to sum into daily values.
    :type values: numpy.ndarray
    :rtype: numpy.ndarray

    """
    new_shape = (-1, 24) + values.shape[1:]
    values = values.reshape(new_shape)
    if minmax:
        return np.sum(values, axis=1), np.max(values, axis=1), np.min(values, axis=1), np.median(values, axis=1), \
               np.average(values, axis=1)
    else:
        return np.sum(values, axis=1)


def bind_hourly_to_daily(values):
    """Create daily values from hourly values by binding the daily values
    in the array.  This is how temperature and ET is currently handled

    :param values: Hourly values already applicable to daily.
    :type values: numpy.ndarray
    :rtype: numpy.ndarray

    """

    new_shape = (-1, 24) + values.shape[1:]
    values = values.reshape(new_shape)

    return np.average(values, axis=1)
