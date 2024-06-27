import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def calculate_rmse(ds1, ds2, total_weight, var='T'):
    # Determine the number of months in ds1
    num_months = ds1[var].shape[0]
    
    # Slice total_weight to match the number of months in ds1
    total_weight_sliced = total_weight[:num_months, :, :]
    
    # Initialize the RMSE array with NaN values
    rmse_per_month = np.full(12, np.nan)
    
    # Compute RMSE for existing months
    squared_diff = (ds1[var] - ds2[var]) ** 2
    weighted_squared_diff = squared_diff * total_weight_sliced
    weighted_sum = weighted_squared_diff.sum(axis=(1, 2))
    total_weight_sum = total_weight_sliced.sum(axis=(1, 2))
    weighted_mean_squared_diff = weighted_sum / total_weight_sum
    rmse_existing_months = np.sqrt(weighted_mean_squared_diff)
    
    # Fill in the RMSE array with the computed values
    rmse_per_month[:num_months] = rmse_existing_months.values
    
    return rmse_per_month

def main(shared_path, hybrid_path_h0):
    # Load the baseline physics simulation datasets
    ds_mmf_ref = xr.open_dataset(shared_path + 'h0/1year/mmf_ref/mmf_ref.eam.h0.0003.nc')
    ds_mmf_a = xr.open_dataset(shared_path + 'h0/1year/mmf_a/mmf_a.eam.h0.0003.nc')
    ds_mmf_b = xr.open_dataset(shared_path + 'h0/1year/mmf_b/mmf_b.eam.h0.0003.nc')
    ds_mmf_c = xr.open_dataset(shared_path + 'h0/1year/mmf_c/mmf_c.eam.h0.0003.nc')
    ds_nn = xr.open_mfdataset(hybrid_path_h0+'*.eam.h0.0003-*.nc', engine='netcdf4')

    # Calculate mass weights
    p_interface = ds_mmf_ref.hyai.values[np.newaxis,:,np.newaxis]*ds_mmf_ref.P0.values + ds_mmf_ref.hybi.values[np.newaxis,:,np.newaxis]*ds_mmf_ref.PS.values[:,np.newaxis,:]
    dp = p_interface[:,1:61,:] - p_interface[:,0:60,:]
    area = ds_mmf_ref.area
    area_weight = area.values[np.newaxis,np.newaxis,:]
    total_weight = dp*area_weight

    # Colors and markers for each group
    colors = ['cyan', 'blue', 'red', 'purple']
    markers = ['o', 's', '^', 'D']

    var1 = 'T'
    var2 = 'Q'
    months = np.arange(1, 13)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Plot for Temperature (T)
    rmse_per_month_nn = calculate_rmse(ds_nn, ds_mmf_ref, total_weight)
    ax1.plot(months, rmse_per_month_nn, label='NN', color='red', linestyle='--', marker='x')

    rmse_per_month_mmf_a = calculate_rmse(ds_mmf_a, ds_mmf_ref, total_weight)
    ax1.plot(months, rmse_per_month_mmf_a, label='MMF', color='black', linestyle='--', marker='x')

    rmse_per_month_mmf_b = calculate_rmse(ds_mmf_b, ds_mmf_ref, total_weight)
    ax1.plot(months, rmse_per_month_mmf_b, color='black', linestyle='--', marker='^')

    rmse_per_month_mmf_c = calculate_rmse(ds_mmf_c, ds_mmf_ref, total_weight)
    ax1.plot(months, rmse_per_month_mmf_c, color='black', linestyle='--', marker='^')

    ax1.set_yscale('log')
    ax1.set_yticks([0.5, 1, 2, 5, 10, 20, 50, 100])
    ax1.set_yticklabels(['0.5', '1', '2', '5', '10', '20', '50', '100'], fontsize=12)
    ax1.set_xlabel('Month', fontsize=15)
    ax1.set_ylabel('online RMSE (K)', fontsize=15)
    ax1.set_title('Temperature root mean squared error (RMSE)', fontsize=15)
    ax1.legend(fontsize=10, loc='upper left')
    ax1.grid(True)
    ax1.set_xticks([0] + list(months))
    ax1.set_xticklabels([0] + list(months), fontsize=12)

    # Plot for Moisture (Q)
    rmse_per_month_nn = calculate_rmse(ds_nn, ds_mmf_ref, total_weight, var='Q')*1e3
    ax2.plot(months, rmse_per_month_nn, label='NN', color='red', linestyle='--', marker='x')

    rmse_per_month_mmf_a = calculate_rmse(ds_mmf_a, ds_mmf_ref, total_weight, var='Q')*1e3
    ax2.plot(months, rmse_per_month_mmf_a, label='MMF', color='black', linestyle='--', marker='x')

    rmse_per_month_mmf_b = calculate_rmse(ds_mmf_b, ds_mmf_ref, total_weight, var='Q')*1e3
    ax2.plot(months, rmse_per_month_mmf_b, color='black', linestyle='--', marker='^')

    rmse_per_month_mmf_c = calculate_rmse(ds_mmf_c, ds_mmf_ref, total_weight, var='Q')*1e3
    ax2.plot(months, rmse_per_month_mmf_c, color='black', linestyle='--', marker='^')

    ax2.set_yscale('log')
    ax2.set_yticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
    ax2.set_yticklabels(['0.1', '0.2', '0.5', '1', '2', '5', '10'], fontsize=12)
    ax2.set_xlabel('Month', fontsize=15)
    ax2.set_ylabel('online RMSE (g/kg)', fontsize=15)
    ax2.set_title('Moisture root mean squared error (RMSE)', fontsize=15)
    # ax2.legend(fontsize=10, loc='upper right')  # Uncomment if needed
    ax2.grid(True)
    ax2.set_xticks([0] + list(months))
    ax2.set_xticklabels([0] + list(months), fontsize=12)

    plt.tight_layout()

    # Ensure the directory exists before saving the plot
    os.makedirs('./figure', exist_ok=True)
    plt.savefig('./figure/downstream_monthly_rmse.jpeg')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot RMSE for Temperature and Moisture")
    parser.add_argument('shared_path', type=str, help='Path to the shared directory')
    parser.add_argument('hybrid_path_h0', type=str, help='Path to the hybrid h0 directory')
    args = parser.parse_args()

    main(args.shared_path, args.hybrid_path_h0)
