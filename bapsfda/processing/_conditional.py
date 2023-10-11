# import numpy as np
# from scipy.signal import find_peaks

# def gather_peaks(flucts, trig_rms, multi, peak_distance):
#     """
#     Finds peaks in a signal above a multiple of the RMS of the signal.

#     Args:
#         flucts (_type_): _description_
#         trig_rms (_type_): _description_
#         multi (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     p_list = []
#     count = 0
#     number_empty = 0
#     empt = []
#     for i, shot in enumerate(flucts):
#         p, _ = find_peaks(shot, height=trig_rms[i] * multi, distance=peak_distance)
#         p_list.append(p)
#         if len(p_list[i]) == 0:
#             number_empty += 1
#             empt.append(i)
#         count += len(p_list[i])
#     print(
#         f"total events: {count}, shots with no events: {number_empty}, empty shots: {empt}"
#     )
#     return p_list


# def slice_event(peaks, signal, left_offset, right_offset):
#     """
#     For each peak in a list, takes a slice of the signal around that peak.

#     Args:
#         peaks (_type_): List of peak indices
#         signal (_type_): Signal to be sliced
#         left_offset (_type_): _description_
#         right_offset (_type_): _description_

#     Returns:
#         _type_: Array of event slices
#     """
#     t_end = signal.shape[-1]
#     b = []
#     for peak_index in peaks:
#         if ((peak_index - left_offset) >= 0) and (
#             (peak_index + right_offset + 1) <= t_end
#         ):
#             b.append(signal[peak_index - left_offset : peak_index + right_offset + 1])
#     if b == []:
#         return np.array([0])
#     return np.array(b)


# def conditional_average(peak_list, flucts, lo, ro):
#     e_list = []
#     for i, shot in enumerate(flucts):
#         b = []
#         empt = 0
#         for m in peak_list[i]:
#             if ((m - lo) >= 0) and ((m + ro + 1) <= t_end):
#                 b.append(shot[m - lo : m + ro + 1])
#         b = np.array(b)
#         try:
#             b.shape[1]
#             a = np.mean(b, axis=0)
#         except IndexError:
#             empt += 1
#         e_list.append(a)
#     return np.array(e_list)
