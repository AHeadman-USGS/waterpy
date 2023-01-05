def main():

    """Main module that runs waterpy.

    This module contains functionality that:
    - Read model configurationo file
    - Read all input files
    - Preprocess input data
        - Calculate the timestep daily fraction
        - Calculate pet if not in timeseries
        - Calculates adjusted precipitation from snowmelt
        - Calculate the twi weighted mean
    - Run Topmodel
    - Post process results
        - Write output *.csv file of results
        - Plot output
        """
import pandas as pd
from pathlib import PurePath
from waterpy import (hydrocalcs,
                     modelconfigfile,
                     parametersfile,
                     timeseriesfile,
                     twifile,
                     plots,
                     report)
from waterpy.topmodel import Topmodel


def waterpy(configfile, options):
    """Read inputs, preprocesse data, run Topmodel, and postprocess
    results, write *.csv outputfiles and make plots.

    :param configfile: The file path to the model config file that
    contains model specifications
    :type configfile: string
    :param options: The options sent from the cli
    :type options: Click.obj
    """
    config_data = modelconfigfile.read(configfile)
    parameters, timeseries, twi, database = read_input_files(config_data)
    preprocessed_data = preprocess(config_data, parameters, timeseries, twi)
    topmodel_data = run_topmodel(config_data, parameters, timeseries, twi, preprocessed_data)
    postprocess(config_data, timeseries, preprocessed_data, topmodel_data)


def read_input_files(config_data):
    """Read input files from model configuration file.

    Returns a tuple of:
        dictionary from parameters file
        pandas.DataFrame from timeseries file
        pandas.DataFrame from twi file

    :param config_data: A ConfigParser object that behaves much like a dictionary.
    :type config_data: ConfigParser
    :return: Tuple of parameters dict, timeseries dataframe, twi dataframe
    :rtype: tuple
    """
    characteristics_basin = parametersfile.read(config_data["Inputs"]["characteristics_basin_file"])
    parameters_land_type = parametersfile.read(config_data["Inputs"]["parameters_land_type_file"])
    parameters = {
        "basin": characteristics_basin,
        "land_type": parameters_land_type,
    }
    timeseries = timeseriesfile.read(config_data["Inputs"]["timeseries_file"])
    twi = twifile.read(config_data["Inputs"]["twi_file"])
    database = config_data["Inputs"]["data_dir"]

    return parameters, timeseries, twi, database


def preprocess(config_data, parameters, timeseries, twi):
    """Preprocess data for topmodel run.

    Calculate timestep daily fraction, usually 1 for daily timesteps
        - 1 day = 86400 seconds
    Calculate pet if pet is not in timeseries dataframe
    Calculate snowmelt and adjusted precipitation from snowmelt routine
        - Snowmelt routine requires temperatures in Fahrenheit.
        - The temperature cutoff from the parameters dict is in Fahrenheit.
        - snowprecip is the adjusted precipitation from snowmelt.
    Calculate the difference between the adjusted precip and pet for Topmodel.
    Calculate the weighted twi mean for Topmodel.

    :param config_data: A ConfigParser object that behaves much like a dictionary.
    :type config_data: ConfigParser
    :param parameters: The parameters for the model.
    :type parameters: Dict
    :param timeseries: A dataframe of all the timeseries data.
    :type timeseries: Pandas.DataFrame
    :param twi: A dataframe of all the twi data.
    :type twi: Pandas.DataFrame
    :return preprocessed_data: A dict of the calculated variables from
                               preprocessing.
    :rtype: dict
    """
    # Calculate the daily timestep as a fraction
    if config_data["Options"].getboolean("option_distribution_record"):
        rain_dist_file = config_data["Options"]["option_dist_file"]
    else:
        rain_dist_file = None

    timestep_daily_fraction = (
        (timeseries.index[1] - timeseries.index[0]).total_seconds() / 86400.0
    )

    # Get pet as a numpy array from the input timeseries if it exists,
    # otherwise calculate it.
    if "pet" in timeseries.columns:
        pet = timeseries["pet"].to_numpy() * timestep_daily_fraction
    else:
        pet = hydrocalcs.pet(
            dates=timeseries.index.to_pydatetime(),
            temperatures=timeseries["temperature"].to_numpy(),
            latitude=parameters["basin"]["latitude"]["value"],
            calib_coeff=parameters["land_type"]["pet_calib_coeff"]["value"],
            method="hamon"
        )
        pet = pet * timestep_daily_fraction

    # If snowmelt option is turned on, then compute snowmelt and the difference
    # between the adjusted precip with pet.
    # Otherwise, just compute the difference between the original precip with
    # pet.
    snowprecip = None
    snowmelt = None
    snowpack = None
    snow_water_equivalence = None
    if config_data["Options"].getboolean("option_snowmelt"):
        # Calculate the adjusted precipitation based on snowmelt
        # Note: snowmelt function needs temperatures in Fahrenheit
        snowprecip, snowmelt, snowpack, snow_water_equivalence = hydrocalcs.snowmelt(
            timeseries["precipitation"].to_numpy(),
            timeseries["temperature"].to_numpy() * (9/5) + 32,
            parameters["land_type"]["snowmelt_temperature_cutoff"]["value"],
            parameters["land_type"]["snowmelt_rate_coeff_with_rain"]["value"],
            parameters["land_type"]["snowmelt_rate_coeff"]["value"],
            timestep_daily_fraction
        )

        # Calculate the difference between the adjusted precip (snowprecip)
        # and pet.
        precip_minus_pet = snowprecip - pet
    else:
        # Calculate the difference between the original precip and pet
        precip_minus_pet = timeseries["precipitation"].to_numpy() - pet
    # Calculate the twi weighted mean
    twi_weighted_mean = hydrocalcs.weighted_mean(values=twi["twi"],
                                                 weights=twi["proportion"]) / parameters["land_type"]["twi_adj"]["value"]

    # Adjust the scaling parameter by the spatial coefficient
    scaling_parameter_adjusted = (
        parameters["basin"]["scaling_parameter"]["value"]
        * parameters["land_type"]["spatial_coeff"]["value"]
    )
    # Define basin area value
    basin_area = parameters["basin"]["basin_area_total"]["value"]

    # Return a dict of calculated data
    preprocessed_data = {
        "timestep_daily_fraction": timestep_daily_fraction,
        "pet": pet,
        "precip_minus_pet": precip_minus_pet,
        "snowprecip": snowprecip,
        "snowmelt": snowmelt,
        "snowpack": snowpack,
        "snow_water_equivalence": snow_water_equivalence,
        "twi_weighted_mean": twi_weighted_mean,
        "scaling_parameter_adjusted": scaling_parameter_adjusted,
        "basin_area": basin_area,
        "rain_file": rain_dist_file
    }

    return preprocessed_data


def run_topmodel(config_data, parameters, timeseries, twi, preprocessed_data):
    """Run Topmodel.

    :param config_data: A ConfigParser object that behaves much like a dictionary.
    :type config_data: ConfigParser
    :param parameters: The parameters for the model.
    :type parameters: Dict
    :param twi: A dataframe of all the twi data.
    :type twi: Pandas.DataFrame
    :param preprocessed_data: A dict of the calculated variables from
                              preprocessing.
    :type: dict
    :param timeseries: A dataframe for timeseries of temperature, precip and observed flow
    :type: timeseries: Pandas.DataFrame
    :return topmodel_data: A dict of relevant data results from Topmodel
    :rtype: dict
    """
    # Initialize Topmodel
    topmodel = Topmodel(
        scaling_parameter=preprocessed_data["scaling_parameter_adjusted"],
        raw_scaling_parameter=parameters["basin"]["scaling_parameter"]["value"],
        saturated_hydraulic_conductivity=(
            parameters["basin"]["saturated_hydraulic_conductivity"]["value"]
        ),
        saturated_hydraulic_conductivity_multiplier=(
            parameters["basin"]["saturated_hydraulic_conductivity_multiplier"]["value"]
        ),
        macropore_fraction=parameters["land_type"]["macropore_fraction"]["value"],
        soil_depth_total=parameters["basin"]["soil_depth_total"]["value"],
        rooting_depth_factor=parameters["land_type"]["rooting_depth_factor"]["value"],
        field_capacity_fraction=parameters["basin"]["field_capacity_fraction"]["value"],
        porosity_fraction=parameters["basin"]["porosity_fraction"]["value"],
        wilting_point_fraction=parameters["basin"]["wilting_point_fraction"]["value"],
        basin_area_total=parameters["basin"]["basin_area_total"]["value"],
        impervious_area_fraction=parameters["basin"]["impervious_area_fraction"]["value"] / 100,
        impervious_curve_number=parameters["land_type"]["impervious_curve_number"]["value"],
        flow_initial=parameters["basin"]["flow_initial"]["value"],
        twi_adj=parameters["land_type"]["twi_adj"]["value"],
        eff_imp=parameters["land_type"]["eff_imp"]["value"],
        et_exp_dorm=parameters["land_type"]["et_exp_dorm"]["value"],
        et_exp_grow=parameters["land_type"]["et_exp_grow"]["value"],
        grow_trigger=parameters["land_type"]["grow_trigger"]["value"],
        riparian_area=parameters["basin"]["rip_area"]["value"],
        twi_values=twi["twi"].to_numpy(),
        twi_saturated_areas=twi["proportion"].to_numpy(),
        twi_mean=preprocessed_data["twi_weighted_mean"],
        precip_available=preprocessed_data["precip_minus_pet"],
        precip=timeseries["precipitation"],
        pet_hamon=preprocessed_data["pet"],
        temperatures=timeseries["temperature"].to_numpy(),
        timestep_daily_fraction=preprocessed_data["timestep_daily_fraction"],
        rain_file=preprocessed_data["rain_file"],
        option_channel_routing=config_data["Options"].getboolean("option_channel_routing"),
        option_karst=config_data["Options"].getboolean("option_karst"),
        option_randomize_daily_to_hourly=config_data["Options"].getboolean("option_randomize_daily_to_hourly"),
        option_min_max=config_data["Options"].getboolean("option_max_min"),
        option_distribution=config_data["Options"].getboolean("option_distribution_record"),
        option_forecast=config_data["Options"].getboolean("forecast")
    )

    # Run Topmodel
    topmodel.run()

    # Return a dict of relevant calculated values
    topmodel_data = {
        "flow_predicted": topmodel.flow_predicted,
        "saturation_deficit_avgs": topmodel.saturation_deficit_avgs,
        "saturation_deficit_locals": topmodel.saturation_deficit_locals,
        "unsaturated_zone_storages": topmodel.unsaturated_zone_storages,
        "infiltration": topmodel.infiltration_array,
        "root_zone_storages": topmodel.root_zone_storages,
        "evaporations": topmodel.evaporations,
        "infiltration_excess": topmodel.infiltration_excess,
        "evaporation_actual": topmodel.evaporation_actual,
        "precip_available": topmodel.precip_available,
        #"q_root": topmodel.q_root,
        #"sub_flow": topmodel.sub_flow,
        "karst_flow": topmodel.karst_flow,
        "imp_flow": topmodel.flow_predicted_impervious,
        "root_zone_avg": topmodel.root_zone_avg,
        "excesses": topmodel.precip_excesses_op,
        "sat_overland_flow": topmodel.pex_flow,
        "return_flow": topmodel.return_flow_totals,
        #"overland_flow": topmodel.overland_flow
    }

    return topmodel_data


def postprocess(config_data, timeseries, preprocessed_data, topmodel_data):
    """Postprocess data for output.

    Output csv files
    Plot timeseries
    """
    # Get output timeseries data
    timeseries = timeseries[365:]
    output_df = get_output_dataframe(timeseries,
                                     preprocessed_data,
                                     topmodel_data)

    # Get output comparison stats
    output_comparison_data = get_comparison_data(output_df, minmax=config_data["Options"].getboolean("option_max_min"))

    # Write output data
    write_output_csv(df=output_df,
                     filename=PurePath(
                         config_data["Outputs"]["output_dir"],
                         config_data["Outputs"]["output_filename"]),
                     minmax=config_data["Options"].getboolean("option_max_min"))

    # Write output data matrices
    if config_data["Options"].getboolean("option_write_output_matrices"):
        write_output_matrices_csv(config_data, timeseries, topmodel_data)

    # Plot output data
    plot_output_data(df=output_df,
                     comparison_data=output_comparison_data,
                     path=config_data["Outputs"]["output_dir"],
                     minmax=config_data["Options"].getboolean("option_max_min"))

    # Write report of output data
    write_output_report(df=output_df,
                        comparison_data=output_comparison_data,
                        filename=PurePath(
                            config_data["Outputs"]["output_dir"],
                            config_data["Outputs"]["output_report"]),
                        minmax=config_data["Options"].getboolean("option_max_min"))


def get_output_dataframe(timeseries, preprocessed_data, topmodel_data):
    """Get the output data of interest.

    Returns a Pandas Dataframe of all output data of interest.
    """
    output_data = {}
    if preprocessed_data["snowprecip"] is not None:
        output_data["snowprecip"] = preprocessed_data["snowprecip"]
        output_data["snowmelt"] = preprocessed_data["snowmelt"]
        output_data["snowpack"] = preprocessed_data["snowpack"]
        output_data["snow_water_equivalence"] = preprocessed_data["snow_water_equivalence"]

    if "pet" not in timeseries.columns:
        output_data["pet"] = preprocessed_data["pet"][365:]
    output_data["aet"] = topmodel_data["evaporation_actual"]
    output_data["precip_minus_pet"] = preprocessed_data["precip_minus_pet"][365:]
    output_data["infiltration"] = topmodel_data["infiltration"]
    output_data["infiltration_excess"] = topmodel_data["infiltration_excess"]
    #output_data["q_root"] = topmodel_data["q_root"]
    #output_data["sub_flow"] = topmodel_data["sub_flow"]
    output_data["karst_flow"] = topmodel_data["karst_flow"]
    output_data["imp_flow"] = topmodel_data["imp_flow"]
    output_data["root_zone_avg"] = topmodel_data["root_zone_avg"]
    if type(topmodel_data["flow_predicted"]) == tuple:
        output_data["flow_predicted_total"] = topmodel_data["flow_predicted"][0]
        output_data["flow_predicted_maximum"] = topmodel_data["flow_predicted"][1]
        output_data["flow_predicted_minimum"] = topmodel_data["flow_predicted"][2]
        output_data["flow_predicted_median"] = topmodel_data["flow_predicted"][3]
        output_data["flow_predicted_average"] = topmodel_data["flow_predicted"][4]
        output_data["discharge_predicted_maximum"] = (
                ((preprocessed_data["basin_area"] * 1000000 * 35.31467) / 3600000) * topmodel_data["flow_predicted"][1]
        )
        output_data["discharge_predicted_minimum"] = (
                ((preprocessed_data["basin_area"] * 1000000 * 35.31467) / 3600000) * topmodel_data["flow_predicted"][2]
        )
        output_data["discharge_predicted_median"] = (
                ((preprocessed_data["basin_area"] * 1000000 * 35.31467) / 3600000) * topmodel_data["flow_predicted"][3]
        )
        output_data["discharge_predicted_average"] = (
                ((preprocessed_data["basin_area"] * 1000000 * 35.31467) / 3600000) * topmodel_data["flow_predicted"][4]
        )
    else:
        output_data["flow_predicted"] = topmodel_data["flow_predicted"]
        m = preprocessed_data["basin_area"]
        n = topmodel_data["flow_predicted"]
        g = m * 1000000 * 35.31467
        d = g / 86400000
        j = d * n
        output_data["discharge_predicted"] = j
    output_data["saturation_deficit_avgs"] = topmodel_data["saturation_deficit_avgs"]
    output_data["sat_overland_flow"] = topmodel_data["sat_overland_flow"]
    output_data["return_flow"] = topmodel_data["return_flow"]
    #output_data["overland_flow"] = topmodel_data["overland_flow"]
    output_df = timeseries.assign(**output_data)

    return output_df


def get_comparison_data(output_df, minmax):
    """Get comparison statistics.

    Return a dictionary of descriptive statistics and if output data contains
    an observed flow, then compute the Nash-Sutcliffe statistic.
    """
    output_comparison_data = {}
    if minmax:
        if "flow_observed" in output_df.columns:
            output_comparison_data["nash_sutcliffe"] = (
                hydrocalcs.nash_sutcliffe(
                    observed=output_df["flow_observed"].to_numpy(),
                    modeled=output_df["flow_predicted_total"].to_numpy())
            )
            output_comparison_data["absolute_error"] = (
                hydrocalcs.absolute_error(
                    observed=output_df["flow_observed"].to_numpy(),
                    modeled=output_df["flow_predicted_total"].to_numpy())
            )
            output_comparison_data["mean_squared_error"] = (
                hydrocalcs.mean_squared_error(
                    observed=output_df["flow_observed"].to_numpy(),
                    modeled=output_df["flow_predicted_total"].to_numpy())
            )
    else:
        if "flow_observed" in output_df.columns:
            output_comparison_data["nash_sutcliffe"] = (
                hydrocalcs.nash_sutcliffe(
                    observed=output_df["flow_observed"].to_numpy(),
                    modeled=output_df["flow_predicted"].to_numpy())
            )
            output_comparison_data["absolute_error"] = (
                hydrocalcs.absolute_error(
                    observed=output_df["flow_observed"].to_numpy(),
                    modeled=output_df["flow_predicted"].to_numpy())
            )
            output_comparison_data["mean_squared_error"] = (
                hydrocalcs.mean_squared_error(
                    observed=output_df["flow_observed"].to_numpy(),
                    modeled=output_df["flow_predicted"].to_numpy())
            )

    return output_comparison_data


def write_output_csv(df, filename, minmax):
    """Write output timeseries to csv file.

    Creating a pandas Dataframe to ease of saving a csv.
    """
    if minmax is False:
        df = df.rename(columns={
            "temperature": "temperature (celsius)",
            "precipitation": "precipitation (mm/day)",
            "pet": "pet (mm/day)",
            "aet": "aet (mm/day)",
            "precip_minus_pet": "precip_minus_pet (mm/day)",
            "flow_observed": "flow_observed (mm/day)",
            "flow_predicted": "flow_predicted (mm/day)",
            "infiltration": "infiltration (mm/day)",
            "discharge_predicted": "discharge_predicted (cfs)",
            "saturation_deficit_avgs": "saturation_deficit_avgs (mm/day)",
            "snowprecip": "snowprecip (mm/day)",
        })
        df.to_csv(filename,
                  float_format="%.16f")

    else:
        df = df.rename(columns={
            "temperature": "temperature (celsius)",
            "precipitation": "precipitation (mm/day)",
            "pet": "pet (mm/day)",
            "aet": "aet (mm/day)",
            "precip_minus_pet": "precip_minus_pet (mm/day)",
            "flow_observed": "flow_observed (mm/day)",
            "flow_predicted_minimum": "flow_predicted minimum (mm/hour)",
            "flow_predicted_maximum": "flow_predicted maximum (mm/hour)",
            "flow_predicted_median": "flow_predicted median (mm/hour)",
            "flow_predicted_average": "flow_predicted average (mm/hour)",
            "flow_predicted_total": "flow_predicted (mm/day)",
            "discharge_predicted_minimum": "discharge_predicted minimum (cfs)",
            "discharge_predicted_maximum": "discharge_predicted maximum (cfs)",
            "discharge_predicted_median": "discharge_predicted median (cfs)",
            "discharge_predicted_average": "discharge_predicted average (cfs)",
            "infiltration": "infiltration (mm/day)",
            "saturation_deficit_avgs": "saturation_deficit_avgs (mm/day)",
            "snowprecip": "snowprecip (mm/day)",
        })
        df.to_csv(filename,
                  float_format="%.16f")


def write_output_matrices_csv(config_data, timeseries, topmodel_data):
    """Write output matrices.

    Matrices are of size: len(timeseries) x len(twi_bins)

    The following are the matrices saved.
         saturation_deficit_locals
         unsaturated_zone_storages
         root_zone_storages
    """
    num_cols = topmodel_data["saturation_deficit_locals"].shape[1]
    header = ["bin_{}".format(i) for i in range(1, num_cols+1)]

    saturation_deficit_locals_df = (
        pd.DataFrame(topmodel_data["saturation_deficit_locals"],
                     index=timeseries.index)
    )

    unsaturated_zone_storages_df = (
        pd.DataFrame(topmodel_data["unsaturated_zone_storages"],
                     index=timeseries.index)
    )

    root_zone_storages_df = (
        pd.DataFrame(topmodel_data["root_zone_storages"],
                     index=timeseries.index)
    )

    evaporations_df = (
        pd.DataFrame(topmodel_data["evaporations"],
                     index=timeseries.index)
    )

    excesses_df = (
        pd.DataFrame(topmodel_data["excesses"],
                     index=timeseries.index)
    )

    saturation_deficit_locals_df.to_csv(
        PurePath(
            config_data["Outputs"]["output_dir"],
            config_data["Outputs"]["output_filename_saturation_deficit_locals"]
        ),
        float_format="%.2f",
        header=header,
    )

    unsaturated_zone_storages_df.to_csv(
        PurePath(
            config_data["Outputs"]["output_dir"],
            config_data["Outputs"]["output_filename_unsaturated_zone_storages"]
        ),
        float_format="%.16f",
        header=header,
    )

    root_zone_storages_df.to_csv(
        PurePath(
            config_data["Outputs"]["output_dir"],
            config_data["Outputs"]["output_filename_root_zone_storages"]
        ),
        float_format="%.16f",
        header=header,
    )

    evaporations_df.to_csv(
        PurePath(
            config_data["Outputs"]["output_dir"],
            config_data["Outputs"]["output_filename_evaporations"]
        ),
        float_format="%.16f",
        header=header,
    )

    excesses_df.to_csv(
        PurePath(
            config_data["Outputs"]["output_dir"],
            config_data["Outputs"]["output_filename_excesses"]
        ),
        float_format="%.16f",
        header=header,
    )
def plot_output_data(df, comparison_data, path, minmax):
    """Plot output timeseries."""
    for key, series in df.iteritems():
        filename = PurePath(path, "{}.png".format(key.split(" ")[0]))
        plots.plot_timeseries(
            dates=df.index.to_pydatetime(),
            values=series.values,
            mean=series.mean(),
            median=series.median(),
            mode=series.mode()[0],
            max=series.max(),
            min=series.min(),
            label="{}".format(key),
            filename=filename)

    if minmax:
        plots.plot_flow_duration_curve(
            values=df["flow_predicted_total"].to_numpy(),
            label="flow_predicted (mm/day)",
            filename=PurePath(path, "flow_duration_curve.png"))

        if "flow_observed" in df.columns:
            plots.plot_timeseries_comparison(
                dates=df.index.to_pydatetime(),
                observed=df["flow_observed"].to_numpy(),
                modeled=df["flow_predicted_total"].to_numpy(),
                absolute_error=comparison_data["absolute_error"],
                nash_sutcliffe=comparison_data["nash_sutcliffe"],
                mean_squared_error=comparison_data["mean_squared_error"],
                label="flow (mm/day)",
                filename=PurePath(path, "flow_observed_vs_flow_predicted.png"))

            plots.plot_flow_duration_curve_comparison(
                observed=df["flow_observed"].to_numpy(),
                modeled=df["flow_predicted_total"].to_numpy(),
                label="flow (mm/day)",
                filename=PurePath(path, "flow_duration_curved_observed_vs_predicted.png"))

    else:
        plots.plot_flow_duration_curve(
            values=df["flow_predicted"].to_numpy(),
            label="flow_predicted (mm/day)",
            filename=PurePath(path, "flow_duration_curve.png"))

        if "flow_observed" in df.columns:
            plots.plot_timeseries_comparison(
                dates=df.index.to_pydatetime(),
                observed=df["flow_observed"].to_numpy(),
                modeled=df["flow_predicted"].to_numpy(),
                absolute_error=comparison_data["absolute_error"],
                nash_sutcliffe=comparison_data["nash_sutcliffe"],
                mean_squared_error=comparison_data["mean_squared_error"],
                label="flow (mm/day)",
                filename=PurePath(path, "flow_observed_vs_flow_predicted.png"))

            plots.plot_flow_duration_curve_comparison(
                observed=df["flow_observed"].to_numpy(),
                modeled=df["flow_predicted"].to_numpy(),
                label="flow (mm/day)",
                filename=PurePath(path, "flow_duration_curved_observed_vs_predicted.png"))


def write_output_report(df, comparison_data, filename, minmax):
    """Write an html web page with interactive plots."""
    plots_html_data = {}
    for key, value in df.iteritems():
        print(key)
        plots_html_data[key] = plots.plot_timeseries_html(
            dates=df.index.to_pydatetime(),
            values=value,
            label="{}".format(key))

    if minmax:
        flow_duration_curve_data = {
            "flow_duration_curve_html": plots.plot_flow_duration_curve_html(
                values=df["flow_predicted_total"].to_numpy(),
                label="flow_predicted (mm/day)")
        }

        if comparison_data:
            comparison_plot_html = plots.plot_timeseries_comparison_html(
                dates=df.index.to_pydatetime(),
                observed=df["flow_observed"].to_numpy(),
                modeled=df["flow_predicted_total"].to_numpy(),
                absolute_error=comparison_data["absolute_error"],
                label="flow (mm/day)")
            comparison_data.update({"comparison_plot_html": comparison_plot_html})

            flow_duration_curve_comparison_hmtl = (
                plots.plot_flow_duration_curve_comparison_html(
                    observed=df["flow_observed"].to_numpy(),
                    modeled=df["flow_predicted_total"].to_numpy(),
                    label="flow (mm/day)")
            )
            flow_duration_curve_data.update(
                {"flow_duration_curve_comparison_html": flow_duration_curve_comparison_hmtl}
            )
    else:
        flow_duration_curve_data = {
            "flow_duration_curve_html": plots.plot_flow_duration_curve_html(
                values=df["flow_predicted"].to_numpy(),
                label="flow_predicted (mm/day)")
        }

        if comparison_data:
            comparison_plot_html = plots.plot_timeseries_comparison_html(
                dates=df.index.to_pydatetime(),
                observed=df["flow_observed"].to_numpy(),
                modeled=df["flow_predicted"].to_numpy(),
                absolute_error=comparison_data["absolute_error"],
                label="flow (mm/day)")
            comparison_data.update({"comparison_plot_html": comparison_plot_html})

            flow_duration_curve_comparison_hmtl = (
                plots.plot_flow_duration_curve_comparison_html(
                    observed=df["flow_observed"].to_numpy(),
                    modeled=df["flow_predicted"].to_numpy(),
                    label="flow (mm/day)")
            )
            flow_duration_curve_data.update(
                {"flow_duration_curve_comparison_html": flow_duration_curve_comparison_hmtl}
            )

    report.save(df=df,
                plots=plots_html_data,
                comparison_data=comparison_data,
                flow_duration_curve_data=flow_duration_curve_data,
                filename=filename)
