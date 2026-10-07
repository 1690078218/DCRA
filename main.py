"""
Get raw BC (i.e. putative barcode) from fastq file

Output:
    1. csv file a the raw BC count
    2. Stats report (stdout):
        a. how many input read in total
        b. how many (# and %) read pass the polyT and adaptor searching process.
        c. # and % with poly T and adaptor find in both ends
        d. # adnd % with no poly T and adpator find
    3. BC rank plot
    4. csv file of the BC whitelist
"""
from collections import defaultdict, Counter
from tqdm import tqdm
import textwrap
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import zipfile
import io
import logging
import gzip
from fast_edit_distance import edit_distance
import logging
from parser import parse_arg
import helper as helper
from config import *
import polyT_adaptor_finder as polyT_adaptor_finder
import read_assignment as read_assignment

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Parse fastq -> polyT_adaptor_finder.Read class
def get_raw_bc_from_reads(reads, min_q=0, kit=None, **kwargs):
    """
    Get putative BC from each reads from a batch of read (can be defined by batch_iterator function)

    Parameters
    ----------
    reads : LIST
        list of reads named tuple with names: id, seq, q_letter
    min_q: INT
        Only count putative bc with minimum value specified
    save_putative_bc: STR
        Output filename for the putative BCs. Will not output anything by default
    Returns
    -------
    1. Counter of high-confidence putative BC
    2. Counter of high-confidence putative BC
    3. pd.DataFrame containing all putative BCs 
    """
    # init
    read_ids = []
    putative_bcs = []
    putative_bc_min_qs = []
    raw_bc = []
    raw_bc_pass = []
    umis = []
    trim_idxs = []
    pre_bc_flankings = []
    post_umi_flankings = []

    for i,r in enumerate(reads):
        
        # create read object 
        read = polyT_adaptor_finder.Read(read_id = r.id, sequence=str(r.seq), 
                    phred_score=r.q_letter, kit=kit, **kwargs)    
        

        read.get_strand_and_raw_bc()
        read_ids.append(read.id)
        putative_bcs.append(read.raw_bc)
        putative_bc_min_qs.append(read.raw_bc_min_q)
        umis.append(read.putative_UMI)
        trim_idxs.append(read.polyT_trimming_idx)
        pre_bc_flankings.append(read.pre_bc_flanking)
        post_umi_flankings.append(read.post_umi_flanking)
        
        if read.raw_bc_min_q and read.raw_bc_min_q >= min_q:     
            raw_bc.append(read.raw_bc)
            
        if read.raw_bc_min_q and read.raw_bc_min_q < min_q:  
            raw_bc_pass.append(100) #tag for low quality putative bc
        else:
            raw_bc_pass.append(read.adaptor_polyT_pass)

    rst_df = pd.DataFrame(
        {'read_id': read_ids,
         'putative_bc': putative_bcs,
         'putative_bc_min_q': putative_bc_min_qs,
        'putative_umi': umis,
        'polyT_end': trim_idxs,
        'pre_bc_flanking': pre_bc_flankings,
        'post_umi_flanking': post_umi_flankings
        }
        )
    return Counter(raw_bc_pass), rst_df

def add_summary(content, args, write_mode='w'):
    with open(args.summary_fn, write_mode) as f:
        f.write(content)
    if not args.minimal_stdout:
        print(content)

def bc_search_qc_report(pass_count, args):
    '''
    Generate report for the putative barcode detection.
    Print stats for
        a. # of input read in total
        b. # and % read pass the polyT and adaptor searching process.
        c. # and % with poly T and adaptor find in both ends
        d. # adnd % with no poly T and adpator find
    Parameters
    ----------
    pass_count : Counter object
        count of different type of error
    min_phred_score :  INT
        min score used to filter out putative BC
    '''
    total_read = sum(pass_count.values())
    
    print_message=textwrap.dedent(
        f'''
        \n----------------------stats of the putative barcodes--------------------------

        Total number of reads: 
            {total_read:,}
        Reads with unambiguous polyT and adapter positions found:            
            {pass_count[0]+ pass_count[100]:,} ({(pass_count[0]+ pass_count[100])/total_read*100:.2f}% of all reads)
            {pass_count[0]:,} in which all bases in the putative BC have Q>={args.minQ}
        Failed Reads: 
            no polyT and adapter positions found: 
                {pass_count[1]:,} ({pass_count[1]/total_read*100:.2f}% of all reads)
            polyT and adapter positions found in both end (fail to determine strand): 
                {pass_count[2]:,} ({pass_count[2]/total_read*100:.2f}% of all reads)
            multiple polyT and adapter found in one end
                {pass_count[10]:,} ({pass_count[10]/total_read*100:.2f}% of all reads)
        -------------------------------------------------------------------------------\n
        ''')
    return print_message

def get_bc_whitelist(raw_bc_count, full_bc_whitelist=None, exp_cells=None, 
                    count_t=None,force_cell_n=None, high_sensitivity_mode=False, 
                    output_empty = True, empty_max_count = np.inf, 
                    out_plot_fn = DEFAULT_KNEE_PLOT_FN, args=None):
    f"""    
    Get a whitelist from all putative cell bc with high-confidence putative bc counts. 
    If the expect number is provided (default), a quantile-based threshold will be 
    calculated to determine the exact cells to be output. Otherwise, a user-specified 
    ount threshold will be used and the cells/Barocdes with counts above the threshold will be output.
    
    If high_sensitivity_mode = True, the high sensitivity (HS) mode is turned on which uses
    more relaxed threshold  

    If in output_empty=True, a list of BCs that are most likely corresponding to 
    empty droplets will also be produced autimatically , which might be useful in 
    downstream analysis.
        Criteria of selecting these BC:
            1. BC in 10x full whitelist, and
            2. At least {out_plot_fn} away from all selected cell-associated BCs in whitelist. and
            3. BC selected as empty droplets should below a certain high-confidence putative barcode count (empty_max_count)

    Args:
        raw_bc_count (Counter): high-confidence putative BC counts for each unique BC
        full_bc_whitelist (str): filename of the 10X whitelist
        exp_cells (int, optional): expected number of cell. Defaults to None.
        count_t (int, optional): count threshold. Defaults to None.

    Raises:
        ValueError: No valid exp_cells or count_t specified

    Returns:
        dict 1: 
            key: barcodes selected
            values: high-confidence putative BC counts
       
        list  (high_sensitivity_mode only):
            barcodes most likely associated to empty droplet
    """

    # Check args:
    if args:
        full_bc_whitelist = args.full_bc_whitelist
        exp_cells = args.expect_cells
        count_t = args.count_threshold
        auto_detect_cells = args.auto_detect_cells
        high_sensitivity_mode = args.high_sensitivity_mode
        out_plot_fn = args.out_plot_fn
        output_empty = args.out_emptydrop_fn
        empty_max_count = args.empty_max_count
        force_cell_n = args.force_cells

        # use the threshold function in config.py
    if high_sensitivity_mode:
        percentile_count_thres = high_sensitivity_threshold_calculation
    else:
        percentile_count_thres = default_count_threshold_calculation

    whole_whitelist = []

    if full_bc_whitelist.endswith('.zip'):
        with zipfile.ZipFile(full_bc_whitelist) as zf:
            # check if there is only 1 file
            assert len(zf.namelist()) == 1

            with io.TextIOWrapper(zf.open(zf.namelist()[0]), 
                                                    encoding="utf-8") as f:
                for line in f:
                    whole_whitelist.append(line.strip())
    else:
        with open(full_bc_whitelist, 'r') as f:
            for line in f:
                whole_whitelist.append(line.strip())
    
    whole_whitelist = set(whole_whitelist)

    raw_bc_count = {k:v for k,v in raw_bc_count.items() if k in whole_whitelist}
    
    # determine real bc based on the count threshold
    # ... existing code ...

    # determine real bc based on the count threshold
    if force_cell_n:
        if force_cell_n > len(raw_bc_count):
            logger.warning(helper.warning_msg(
                f"force_cells ({force_cell_n}) is larger than the number of unique barcodes found in the data ({len(raw_bc_count)})."))
            count_t = 0
        # count threshold is the minimum count of the top N cells
        else:
            count_t = sorted(list(raw_bc_count.values()))[-force_cell_n]
        rank_threshold = force_cell_n

    if count_t is not None:
        # Calculate rank_threshold if not already set
        if 'rank_threshold' not in locals():
            sorted_counts = sorted(raw_bc_count.values(), reverse=True)
            rank_threshold = sum(1 for c in sorted_counts if c >= count_t)

        knee_plot(list(raw_bc_count.values()), count_t, out_plot_fn, raw_bc_count,
                  auto_detected_rank=rank_threshold, show_annotations=True)
        cells_bc = {k: v for k, v in raw_bc_count.items() if v >= count_t}
        if not output_empty:
            return cells_bc, []
        else:
            # create a confident list of empty drops in high sensitivity mode
            logger.info("Creating emtpy droplets barocde list...")
            ept_bc = []
            ept_bc_max_count = min(cells_bc.values())
            ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
            ept_bc_candidate = \
                [k for k, v in raw_bc_count.items() if v < ept_bc_max_count]
            for k in ept_bc_candidate:
                if min([edit_distance(k, x, max_ed=DEFAULT_EMPTY_DROP_MIN_ED) for x in
                        cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
                    ept_bc.append(k)

                # we don't need too much BC in this list
                if len(ept_bc) > DEFAULT_EMPTY_DROP_NUM:
                    break
            return cells_bc, ept_bc


    elif auto_detect_cells:

        # 使用新的 knee point 检测算法

        import pandas as pd

        import numpy as np

        # 准备数据

        sorted_counts = sorted(raw_bc_count.values(), reverse=True)

        rank = np.arange(1, len(sorted_counts) + 1).astype(float)

        count = np.array(sorted_counts).astype(float)

        # 防止 log(0)

        eps = 1e-10

        count = np.clip(count, eps, None)

        # log-log 变换

        log_rank = np.log10(rank)

        log_count = np.log10(count)

        # 计算局部导数

        d = (log_count[1:] - log_count[:-1]) / (log_rank[1:] - log_rank[:-1])

        # 滑动窗口平滑

        window_radius = 5

        d_smooth = np.full_like(d, np.nan)

        for i in range(window_radius, len(d) - window_radius):
            d_smooth[i] = np.mean(d[i - window_radius: i + window_radius + 1])

        # 测试多个导数阈值

        thresholds = list(range(-1, -51, -1))

        threshold_results = []

        for threshold_derivative in thresholds:

            found_threshold = False

            knee_idx = None

            for i in range(len(d_smooth)):

                if not np.isnan(d_smooth[i]) and d_smooth[i] < threshold_derivative:
                    knee_idx = i

                    found_threshold = True

                    break

            if not found_threshold:

                if not np.all(np.isnan(d_smooth)):

                    knee_idx = np.nanargmin(d_smooth)

                    found_threshold = True

                else:

                    knee_idx = -1

            if knee_idx >= 0 and knee_idx < len(d_smooth):

                threshold_rank = rank[knee_idx + 1]

                threshold_count = count[knee_idx + 1]

                smoothed_derivative = d_smooth[knee_idx]

            else:

                threshold_rank = np.nan

                threshold_count = np.nan

                smoothed_derivative = np.nan

            threshold_results.append({

                'Derivative_Threshold': threshold_derivative,

                'Threshold_Rank': threshold_rank,

                'Threshold_Count': threshold_count,

                'Actual_Smoothed_Derivative': smoothed_derivative,

                'Found': found_threshold

            })

        # 计算变化率

        threshold_results_df = pd.DataFrame(threshold_results)

        threshold_results_df['Rank_Change_Rate'] = np.nan

        threshold_results_df['Count_Change_Rate'] = np.nan

        threshold_results_df['Derivative_Change_Rate'] = np.nan

        threshold_results_df['Rank_Absolute_Change'] = np.nan

        threshold_results_df['Count_Absolute_Change'] = np.nan

        threshold_results_df['Relative_Change_Rate'] = np.nan

        for i in range(1, len(threshold_results_df)):

            prev_row = threshold_results_df.iloc[i - 1]

            curr_row = threshold_results_df.iloc[i]

            if prev_row['Found'] and curr_row['Found']:

                rank_change = curr_row['Threshold_Rank'] - prev_row['Threshold_Rank']

                count_change = curr_row['Threshold_Count'] - prev_row['Threshold_Count']

                derivative_change = curr_row['Actual_Smoothed_Derivative'] - prev_row['Actual_Smoothed_Derivative']

                if prev_row['Threshold_Rank'] != 0:

                    rank_change_rate = rank_change / prev_row['Threshold_Rank'] * 100

                else:

                    rank_change_rate = np.nan

                if prev_row['Threshold_Count'] != 0:

                    count_change_rate = count_change / prev_row['Threshold_Count'] * 100

                else:

                    count_change_rate = np.nan

                if prev_row['Actual_Smoothed_Derivative'] != 0:

                    derivative_change_rate = derivative_change / prev_row['Actual_Smoothed_Derivative'] * 100

                else:

                    derivative_change_rate = np.nan

                if not np.isnan(rank_change_rate) and rank_change_rate != 0:

                    relative_change_rate = count_change_rate / rank_change_rate

                else:

                    relative_change_rate = np.nan

                threshold_results_df.loc[i, 'Rank_Change_Rate'] = rank_change_rate

                threshold_results_df.loc[i, 'Count_Change_Rate'] = count_change_rate

                threshold_results_df.loc[i, 'Derivative_Change_Rate'] = derivative_change_rate

                threshold_results_df.loc[i, 'Rank_Absolute_Change'] = rank_change

                threshold_results_df.loc[i, 'Count_Absolute_Change'] = count_change

                threshold_results_df.loc[i, 'Relative_Change_Rate'] = relative_change_rate

        # 计算导数值下降 1 所需的 Rank 变化

        required_rank_changes = []

        for i in range(1, len(threshold_results_df)):

            prev_row = threshold_results_df.iloc[i - 1]

            curr_row = threshold_results_df.iloc[i]

            if prev_row['Found'] and curr_row['Found'] and pd.notna(prev_row['Threshold_Rank']) and pd.notna(
                    curr_row['Threshold_Rank']):
                rank_change = curr_row['Threshold_Rank'] - prev_row['Threshold_Rank']

                derivative_change = curr_row['Derivative_Threshold'] - prev_row['Derivative_Threshold']

                required_rank_changes.append({

                    'Derivative_Range': f"{prev_row['Derivative_Threshold']} to {curr_row['Derivative_Threshold']}",

                    'Derivative_Change': derivative_change,

                    'Rank_Change': rank_change,

                    'Rank_Change_Per_Derivative_Unit': rank_change,

                    'Start_Rank': prev_row['Threshold_Rank'],

                    'End_Rank': curr_row['Threshold_Rank'],

                    'Start_Count': prev_row['Threshold_Count'],

                    'End_Count': curr_row['Threshold_Count']

                })

        # 检测最大值之前的区域

        selected_threshold = None

        closest_row = None

        if required_rank_changes:

            rank_change_df = pd.DataFrame(required_rank_changes)

            rank_change_series = rank_change_df['Rank_Change_Per_Derivative_Unit'].values

            if len(rank_change_series) > 0:

                max_idx = np.argmax(rank_change_series)

                max_value = rank_change_series[max_idx]

                if max_idx == 0:
                    logger.info("Maximum value is at the first position → search using the full range")
                    pre_max_df = rank_change_df.copy()  # 使用完整数据
                else:
                    logger.info(f"The maximum value is at position  {max_idx}  → only analyze the region before the maximum value")
                    pre_max_df = rank_change_df.iloc[:max_idx].copy()

                    # 检测低值点
                low_threshold = np.max(rank_change_series) * 0.05
                min_continuous = 2
                # 块长上限：仅当最大值位于首区间（max_idx == 0，全范围搜索）时生效（Algorithm S1 第 11 行）
                max_continuous = 50

                pre_max_df.loc[:, 'Is_Low'] = pre_max_df['Rank_Change_Per_Derivative_Unit'] <= low_threshold

                # 查找连续低值群
                continuous_groups = []
                current_start = None
                current_length = 0

                for i in range(len(pre_max_df)):
                    if pre_max_df.iloc[i]['Is_Low']:
                        if current_start is None:
                            current_start = i
                        current_length += 1
                    else:
                        if current_length >= min_continuous:
                            if max_idx == 0:
                                effective_length = min(current_length, max_continuous)
                            else:
                                effective_length = current_length
                            group_rows = pre_max_df.iloc[current_start:current_start + effective_length]
                            midpoints = []
                            for _, row in group_rows.iterrows():
                                s, e = map(float, row['Derivative_Range'].split(' to '))
                                midpoints.append((s + e) / 2)
                            avg_deriv = np.mean(midpoints)
                            continuous_groups.append({
                                'start_idx': current_start,
                                'length': effective_length,
                                'avg_derivative': avg_deriv,
                                'first_range': group_rows.iloc[0]['Derivative_Range'],
                                'last_range': group_rows.iloc[-1]['Derivative_Range'],
                                'strategy': 'continuous_group'
                            })
                        current_start = None
                        current_length = 0

                # 处理最后一个群
                if current_length >= min_continuous:
                    if max_idx == 0:
                        effective_length = min(current_length, max_continuous)
                    else:
                        effective_length = current_length
                    group_rows = pre_max_df.iloc[current_start:current_start + effective_length]
                    midpoints = []
                    for _, row in group_rows.iterrows():
                        s, e = map(float, row['Derivative_Range'].split(' to '))
                        midpoints.append((s + e) / 2)
                    avg_deriv = np.mean(midpoints)
                    continuous_groups.append({
                        'start_idx': current_start,
                        'length': effective_length,
                        'avg_derivative': avg_deriv,
                        'first_range': group_rows.iloc[0]['Derivative_Range'],
                        'last_range': group_rows.iloc[-1]['Derivative_Range'],
                        'strategy': 'continuous_group'
                    })

                # 策略选择
                selected_group = None
                strategy_used = ""

                if continuous_groups:
                    first_group = continuous_groups[0]
                    strategy_used = "continuous_group"
                    selected_group = first_group
                    selected_threshold = round(first_group['avg_derivative'], 1)
                else:
                    if len(pre_max_df) > 0:
                        min_val = np.min(pre_max_df['Rank_Change_Per_Derivative_Unit'])
                        min_indices = np.where(pre_max_df['Rank_Change_Per_Derivative_Unit'] == min_val)[0]
                        if len(min_indices) > 0:
                            first_min_idx = min_indices[0]
                            min_row = pre_max_df.iloc[first_min_idx]
                            derivative_range_str = min_row['Derivative_Range']
                            s, e = map(float, derivative_range_str.split(' to '))
                            if s < e:
                                selected_threshold = s
                            else:
                                selected_threshold = e
                            strategy_used = "first_min_use_earlier"
                            selected_group = {
                                'start_idx': first_min_idx,
                                'length': 1,
                                'avg_derivative': selected_threshold,
                                'selected_threshold': selected_threshold,
                                'first_range': derivative_range_str,
                                'last_range': derivative_range_str,
                                'strategy': 'first_min_use_earlier'
                            }

                # 查找最接近的阈值
                if selected_threshold is not None:
                    exact_match_mask = threshold_results_df['Derivative_Threshold'] == selected_threshold
                    if exact_match_mask.any():
                        closest_row = threshold_results_df[exact_match_mask].iloc[0]
                    else:
                        threshold_results_df['abs_diff'] = np.abs(
                            threshold_results_df['Derivative_Threshold'] - selected_threshold)
                        closest_row = threshold_results_df.loc[threshold_results_df['abs_diff'].idxmin()]

                # DCRA 最终阈值：连续低值块内最后一个区间的 End_Rank（论文 Algorithm S1 第 12 行）
                if strategy_used == "continuous_group":
                    last_idx = selected_group['start_idx'] + selected_group['length'] - 1
                    rank_threshold = int(rank_change_df.iloc[last_idx]['End_Rank'])
                elif closest_row is not None:
                    rank_threshold = int(closest_row['Threshold_Rank'])
                else:
                    rank_threshold = None

                # 如果没有找到合适的阈值，使用默认方法
                if rank_threshold is None or rank_threshold <= 0:
                    t = percentile_count_thres(list(raw_bc_count.values()), exp_cells)
                    knee_plot(list(raw_bc_count.values()), t, out_plot_fn, raw_bc_count,
                              auto_detected_rank=None, show_annotations=False)
                    cells_bc = {k: v for k, v in raw_bc_count.items() if v > t}
                else:
                    # 使用 DCRA 秩阈值：按计数降序保留秩 <= rank_threshold 的条形码（与论文口径一致）
                    sorted_items = sorted(raw_bc_count.items(), key=lambda kv: kv[1], reverse=True)
                    rank_threshold = min(rank_threshold, len(sorted_items))
                    cells_bc = dict(sorted_items[:rank_threshold])
                    t = float(sorted_items[rank_threshold - 1][1])
                    logger.info(f"Knee point detected using derivative change-rate analysis (DCRA):")
                    logger.info(f"  Strategy: {strategy_used}")
                    if closest_row is not None:
                        logger.info(f"  Derivative threshold: {closest_row['Derivative_Threshold']}")
                    logger.info(f"  Rank threshold: {rank_threshold:,}")
                    logger.info(f"  Count threshold: {t:.1f}")
                    knee_plot(list(raw_bc_count.values()), t, out_plot_fn, raw_bc_count,
                              auto_detected_rank=rank_threshold, show_annotations=True)

                if not output_empty:
                    return cells_bc, []
                else:
                    logger.info("Creating emtpy droplets barocde list...")
                    ept_bc = []
                    ept_bc_max_count = min(cells_bc.values())
                    ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
                    ept_bc_candidate = \
                        [k for k, v in raw_bc_count.items() if v < ept_bc_max_count]
                    for k in ept_bc_candidate:
                        if min([edit_distance(k, x, max_ed=DEFAULT_EMPTY_DROP_MIN_ED) for x in
                                cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
                            ept_bc.append(k)
                        # we don't need too much BC in this list
                        if len(ept_bc) > DEFAULT_EMPTY_DROP_NUM:
                            break
                    return cells_bc, ept_bc
            else:
                raise ValueError('Invalid value of count_t and/or auto_detect_cells.')

def knee_plot(counts, threshold=None, out_fn='knee_plot.png', barcode_counts_dict=None,
              auto_detected_rank=None, show_annotations=False):
    """
    Plot knee plot using the high-confidence putative BC counts and output sorted barcode triplets

    Args:
        counts (list): high-confidence putative BC counts
        threshold (int, optional): a line to show the count threshold. Defaults to None.
        out_fn (str): output filename for the knee plot. Defaults to 'knee_plot.png'
        barcode_counts_dict (dict, optional): dictionary mapping barcodes to their counts. If provided,
                                            a CSV file with sorted barcodes will be generated.
        auto_detected_rank (int, optional): automatically detected rank threshold. Defaults to None.
        show_annotations (bool, optional): whether to show rank and count annotations on the plot. Defaults to False.
    """
    counts = sorted(counts)[::-1]

    # Create the figure
    plt.figure(figsize=(10, 8))
    plt.loglog(counts, marker='o', linestyle="", alpha=0.6, markersize=4, color='steelblue')

    # Add threshold line
    if threshold is not None:
        plt.axhline(y=threshold, color='red', linestyle='--', linewidth=2, label=f'Count threshold: {threshold:.0f}')

        # Add vertical line at the rank position if auto_detected_rank is provided
        if auto_detected_rank is not None:
            plt.axvline(x=auto_detected_rank, color='green', linestyle='--', linewidth=2,
                        label=f'Rank threshold: {auto_detected_rank:,}')

            # Mark the knee point
            if auto_detected_rank <= len(counts):
                knee_count = counts[auto_detected_rank - 1]
                plt.plot(auto_detected_rank, knee_count, 'ro', markersize=12, zorder=5)
                plt.plot(auto_detected_rank, knee_count, 'r+', markeredgewidth=3, markersize=15, zorder=6)

    # Set labels and title
    plt.xlabel('Barcode Rank', fontsize=12, fontweight='bold')
    plt.ylabel('Read Counts', fontsize=12, fontweight='bold')
    plt.title('Barcode Rank Plot (Knee Point Detection)', fontsize=14, fontweight='bold')

    # Add annotations if requested
    if show_annotations and threshold is not None:
        # Find the rank corresponding to the threshold
        if auto_detected_rank is not None:
            rank_at_threshold = auto_detected_rank
            count_at_threshold = threshold

            # Add text box with detailed information
            info_text = f'Auto-Detected Knee Point\n'
            info_text += f'Rank: {rank_at_threshold:,}\n'
            info_text += f'Count: {count_at_threshold:.0f}'

            plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes,
                     fontsize=11, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, edgecolor='orange'))

            # Annotate the knee point
            if auto_detected_rank <= len(counts):
                plt.annotate(f'({rank_at_threshold:,}, {count_at_threshold:.0f})',
                             xy=(auto_detected_rank, knee_count),
                             xytext=(auto_detected_rank * 1.5, knee_count * 1.5),
                             fontsize=10, fontweight='bold',
                             arrowprops=dict(arrowstyle='->', color='red', lw=2, connectionstyle='arc3,rad=0.3'),
                             color='red')

    # Add grid
    plt.grid(True, alpha=0.3, linestyle=':', linewidth=0.7)

    # Add legend
    plt.legend(loc='upper right', fontsize=10, framealpha=0.9)

    # Save the figure
    plt.tight_layout()
    plt.savefig(out_fn, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Knee plot saved to {out_fn}")

    # Output sorted barcodes triplets file if barcode_counts_dict is provided
    if barcode_counts_dict is not None:
        # Create the triplets file name based on the plot file name
        triplets_fn = out_fn.replace('.png', '_triplets.csv')

        # Sort the barcodes by their counts in descending order
        sorted_barcodes = sorted(barcode_counts_dict.items(), key=lambda x: x[1], reverse=True)

        # Write to CSV file
        import csv
        with open(triplets_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write header
            writer.writerow(['Rank', 'Barcode', 'Count'])
            # Write sorted barcodes with their ranks and counts
            for rank, (barcode, count) in enumerate(sorted_barcodes, start=1):
                writer.writerow([rank, barcode, count])

        print(f"Sorted barcode triplets saved to {triplets_fn}")
        print(f"Total barcodes in triplets: {len(sorted_barcodes):,}")

        # Print additional statistics
        if threshold is not None and auto_detected_rank is not None:
            print(f"\nStatistics:")
            print(f"  Barcodes above threshold: {auto_detected_rank:,}")
            print(f"  Minimum count for selected barcodes: {threshold:.0f}")

def read_batch_generator(fastq_fns, batch_size):
    """Generator of barches of reads from list of fastq files

    Args:aipfile
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq = helper.fastq_parser(handle)
                read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            with open(fn, "r") as handle:
                fastq = helper.fastq_parser(handle)
                read_batch = helper.batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch

def print_logo(args):
    if not args.minimal_stdout:
        print(textwrap.dedent(
            f'''\n\nWelcome to
                {DCRA_LOGO}
        '''))

def main(args_string=None):
    # Start running: Welcome logo
    args = parse_arg(args_string)
    print_logo(args)

    # TMP: print all the arguments in args
    #for arg in vars(args):
    #    print(f"{arg}: {getattr(args, arg)}")
    
    ######################
    ###### Getting putative barcodes
    ######################
    if args.do_bc_search:
        logger.info(f'Getting putative barcodes from {len(args.fastq_fns)} FASTQ files...')
        read_batchs = read_batch_generator(args.fastq_fns, batch_size=args.batch_size)

            
        rst_futures = helper.multiprocessing_submit(get_raw_bc_from_reads,
                                                read_batchs, n_process=args.threads, 
                                                min_q=args.minQ, kit=args.kit_version,
                                                umi_len=args.umi_len)
    

        raw_bc_pass_count = defaultdict(int)

        for idx, f in enumerate(rst_futures):
            count_pass, rst_df = f.result() #write rst_df out
            for k,v in count_pass.items():
                raw_bc_pass_count[k] += v
            if idx == 0:
                rst_df.to_csv(args.out_raw_bc_fn, index=False)
            else:
                rst_df.to_csv(args.out_raw_bc_fn, mode='a', index=False, header=False)
        
        helper.green_msg(f'Putative barcode table saved in {args.out_raw_bc_fn}', printit=True)
        
        # ----------------------stats of the putative barcodes--------------------------

        add_summary(bc_search_qc_report(raw_bc_pass_count, args), 
                    args=args, write_mode='w')

    ######################
    ###### Whitelisting
    ######################
    if args.do_whitelisting:
        # get bc count dict (filtered by minQ)
        dfs = pd.read_csv(args.out_raw_bc_fn, chunksize=1_000_000)
        raw_bc_count = Counter()
        for df in tqdm(dfs, desc = 'Counting high-quality putative BC', unit='M reads'):
            raw_bc_count += Counter(df[
                df.putative_bc_min_q >=args.minQ].putative_bc.value_counts().to_dict())

        try:
            bc_whitelist, ept_bc = get_bc_whitelist(raw_bc_count, args=args)

            # output
            # Sort barcodes by count (descending) to match the order in triplets file
            sorted_whitelist = sorted(bc_whitelist.items(), key=lambda x: x[1], reverse=True)

            with open(args.out_whitelist_fn, 'w') as f:
                for barcode, count in sorted_whitelist:
                    f.write(barcode + '\n')
            with open(args.out_emptydrop_fn, 'w') as f:
                for k in ept_bc:
                    f.write(k + '\n')
                helper.green_msg(f'Empty droplet barcode list saved as `{args.out_emptydrop_fn}`.', printit=True)
                # write to summary
            add_summary(f'\nIdentified # of cells: {len(bc_whitelist)}\n', args, write_mode='a')

        except Exception as e:
            logger.exception(e)
            helper.err_msg(
                "Error: Failed to get whitelist. Please check the input files and settings." , printit=True
                )
        else:
            helper.green_msg(f'Whitelist saved as `{args.out_whitelist_fn}`!', printit=True)

        ######################
        ###### Demultiplexing
        ######################
    if args.do_demultiplexing:
        logger.info("Assigning reads to whitelist.\n")
        # write to fastq
        demul_count_tot, count_tot = read_assignment.assign_read(args=args)
        # write to summary
        add_summary(f'\nTotal reads: {count_tot}'
                    f'\nTotal reads in cells: {demul_count_tot} ({demul_count_tot / count_tot * 100:.2f}%)', args,
                    write_mode='a')

if __name__ == '__main__':
    main()
