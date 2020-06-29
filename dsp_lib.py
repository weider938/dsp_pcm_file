import numpy as np
import time
import math
import struct

"""
notes:
float16_value = np.float16(231.12)
float16_bytearray = float16_value.tobytes()
bytes_to_float_value = np.frombuffer(float16_bytearray, dtype=np.float16)[0]
bits_string = bin(np.float16(float16_value).view('H'))[2:].zfill(16)

bytearray_sample = bytearray(struct.pack("f", value_float32))
"""


def change_ptr_in_file(f_name=None):
    with open(f_name + "_zpt.txt", "w") as f:
        f.write("")
    with open(f_name, "r") as f:
        lines_file = f.read()
        lines = lines_file.split("\n")
        for line in lines:
            try:
                list = line.split(".")
                new_line = list[0] + "," + list[1]
            except:
                pass
            with open(f_name + "_zpt.txt", "a") as f2:
                f2.write(new_line)


def open_bin_file(f_name="333.pcm"):
    start_time_read_file = time.time()
    data = np.fromfile(f_name, dtype=np.short)  # np.short == ">H"
    #  data = np.memmap(f_name, dtype=np.short, mode="r")  # np.short == ">H"  ### 2 вариант считывания
    count_chunks = len(data)
    print("[def] File read. Samples total: ", count_chunks)
    print("--- %s seconds ---" % (time.time() - start_time_read_file))
    return data, count_chunks


def fft_numpy(signal, window, output):
    print("\n[def] Calculating spectrum values using the Numpy library")
    start_time = time.time()
    string_to_file = ""
    signal_array_cplx = []
    for i in range(0, window):  # rows
        sign_cplx = complex(signal[i * 2], signal[i * 2 + 1])
        signal_array_cplx.append(sign_cplx)
    Spectr = np.fft.fft(signal_array_cplx)
    coefficients_abs = []
    for i in Spectr:
        value_to_file = float("{0:.3f}".format(abs(i / window)))
        coefficients_abs.append(value_to_file)
        string_to_file += (str(value_to_file) + "\n")
    with open(output, 'w', encoding="utf-8") as file:
        file.write(string_to_file)
    print(" --- %s seconds ---" % (time.time() - start_time))

    return coefficients_abs


def slow_ft(signal=None,
            window=None,
            output=None):
    print("\n[def] Saving DFT coefficients in" + output)
    start_time = time.time()
    c = 0
    string_to_file = ""
    coeff_spectrum_abs = []

    for j in range(0, window):
        x = complex(0)
        for i in range(0, window):  # rows
            sign_cplx = complex(signal[i * 2], signal[i * 2 + 1])
            W = complex(math.cos(float(np.pi * 2 * i * j / window)),
                        float(-math.sin(np.pi * 2 * i * j / window)))
            x = x + (sign_cplx * W)
        c += 1
        z = complex((x.real / window), (x.imag / window))
        value_to_file = float("{0:.3f}".format(abs(z)))
        string_to_file += (str(value_to_file) + "\n")
        coeff_spectrum_abs.append(value_to_file)

    with open(output, 'w', encoding="utf-8") as file:
        file.write(string_to_file)
    print("--- %s seconds ---" % (time.time() - start_time))
    print("[def] Recording " + output + " completed.\n")
    return coeff_spectrum_abs


def ReturnInvInt(my_int, step=13):
    if my_int == 0:
        return 0
    else:
        cel_part = my_int
        return_int = 0
        ost = -1
        ost_arr = []
        while cel_part != 1:
            ost = cel_part % 2
            cel_part = cel_part // 2
            ost_arr.append(ost)
        ost_arr.append(1)
        if len(ost_arr) < step:
            for i in range(0, step - ((len(ost_arr)))):
                ost_arr.append(0)
        cnt = 1
        for y in ost_arr:
            if y == 1:
                greed_2 = len(ost_arr) - cnt
                return_int = return_int + 2 ** greed_2
            cnt = cnt + 1
        return return_int


def GetArrayFFTIndexes(win, step):
    array_to_ret = []
    for i in range(0, win):
        index = ReturnInvInt(my_int=i, step=step)
        array_to_ret.append(index)
    return array_to_ret


def fast_ft(signal=None,
            window=None,
            output=None,
            print_events=True):
    if print_events:
        print("\n[def] Сохранение коэффициентов ДБПФ в " + output)
    start_time = time.time()
    Spectrum_FFT = []
    for i in range(0, window):
        Spectrum_FFT.append(float(0))

    step = 0
    if window == 128:
        step = 7
    elif window == 256:
        step = 8
    elif window == 512:
        step = 9
    elif window == 1024:
        step = 10
    elif window == 2048:
        step = 11
    elif window == 4096:
        step = 12
    elif window == 8192:
        step = 13
    elif window == 16384:
        step = 14
    elif window == 32768:
        step = 15
    elif window == 65536:
        step = 16
    elif window == 131072:
        step = 17
    elif window == 262144:
        step = 18
    else:
        step = -1
        if print_events:
            print(" Window value must be a multiple of degree 2")
    if step != -1:
        mas_index = GetArrayFFTIndexes(win=window, step=step)
        for i in range(0, step):
            size_block = 2 * pow(2, i)
            count_block = int(window / size_block)
            for m in range(0, int(count_block)):
                for k in range(0, int(size_block / 2)):
                    if i == 0:
                        index_0 = int(mas_index[m * size_block])
                        index_1 = int(mas_index[m * size_block + int(size_block / 2)])
                        S_0 = complex((signal[index_0]).real, (signal[index_0]).imag)
                        S_1 = complex((signal[index_1]).real, (signal[index_1]).imag)

                        W = complex(math.cos((3.141592 * 2 / size_block) * k),
                                    math.sin((3.141592 * 2 / size_block) * k))

                        index_0 = int(m * size_block + k)
                        index_1 = int(m * size_block + (size_block / 2) + k)
                        Spectrum_FFT[index_0] = S_0 + W * S_1
                        Spectrum_FFT[index_1] = S_0 - W * S_1
                    else:
                        index_0 = int(m * size_block + k)
                        index_1 = int(m * size_block + (size_block / 2) + k)

                        S_0 = Spectrum_FFT[index_0]
                        S_1 = Spectrum_FFT[index_1]
                        W = complex(math.cos((3.141592 * 2 / size_block) * k),
                                    math.sin((3.141592 * 2 / size_block) * k))

                        Spectrum_FFT[index_0] = S_0 + W * S_1
                        Spectrum_FFT[index_1] = S_0 - W * S_1
        string_to_file = ""
        list_coef_abs = []
        list_coef = []
        for coef in Spectrum_FFT:
            z = complex((coef.real / window), (coef.imag / window))
            list_coef.append(z)
            value_to_file = float("{0:.3f}".format(abs(z)))
            list_coef_abs.append(value_to_file)
        for it in list_coef_abs[::-1]:
            iter = str(it)
            string_to_file += (iter + "\n")
        if output != None:
            with open(output, 'w', encoding="utf-8") as file:
                file.write(string_to_file)
        if print_events:
            print("--- %s seconds ---" % (time.time() - start_time))
            print("")
        return list_coef[::-1], list_coef_abs[::-1]


def fast_ft_reverse(coefficients=None,
                    window=None,
                    print_events=True):
    if print_events:
        print("\n [def] Calculating the values of the time coefficients ODPF")
    start_time = time.time()

    Samples_RFFT = []
    for i in range(0, window):
        Samples_RFFT.append(float(0))

    step = 0
    if window == 128:
        step = 7
    elif window == 256:
        step = 8
    elif window == 512:
        step = 9
    elif window == 1024:
        step = 10
    elif window == 2048:
        step = 11
    elif window == 4096:
        step = 12
    elif window == 8192:
        step = 13
    elif window == 16384:
        step = 14
    elif window == 32768:
        step = 15
    elif window == 65536:
        step = 16
    elif window == 131072:
        step = 17
    elif window == 262144:
        step = 18
    else:
        if print_events:
            print(" ! Window value must be a multiple of degree 2")

    mas_index = GetArrayFFTIndexes(win=window, step=step)

    for i in range(0, step):
        size_block = 2 * pow(2, i)
        count_block = int(window / size_block)
        for m in range(0, int(count_block)):
            for k in range(0, int(size_block / 2)):
                if i == 0:
                    index_0 = int(mas_index[m * size_block])
                    index_1 = int(mas_index[m * size_block + int(size_block / 2)])
                    S_0 = complex((coefficients[index_0]).real, (coefficients[index_0]).imag)
                    S_1 = complex((coefficients[index_1]).real, (coefficients[index_1]).imag)

                    W = complex(math.cos((3.141592 * 2 / size_block) * k), math.sin((3.141592 * 2 / size_block) * k))

                    index_0 = int(m * size_block + k)
                    index_1 = int(m * size_block + (size_block / 2) + k)
                    Samples_RFFT[index_0] = S_0 + W * S_1
                    Samples_RFFT[index_1] = S_0 - W * S_1
                else:
                    index_0 = int(m * size_block + k)
                    index_1 = int(m * size_block + (size_block / 2) + k)

                    S_0 = Samples_RFFT[index_0]
                    S_1 = Samples_RFFT[index_1]
                    W = complex(math.cos((3.141592 * 2 / size_block) * k), math.sin((3.141592 * 2 / size_block) * k))

                    Samples_RFFT[index_0] = S_0 + W * S_1
                    Samples_RFFT[index_1] = S_0 - W * S_1

    Samples_RFFT_abs = []

    for sps in Samples_RFFT:
        z = complex((sps.real), (sps.imag))
        value = float("{0:.3f}".format(abs(z)))
        Samples_RFFT_abs.append(value)

    if print_events:
        print("--- %s seconds ---" % (time.time() - start_time))
        print("")
    return Samples_RFFT, Samples_RFFT_abs


def return_samples(signal=None,
                   window=None):
    sign_cplx = []
    for i in range(0, window):
        sign_cplx.append(complex(signal[i * 2], signal[i * 2 + 1]))

    sign_cplx_abs = []
    for h in range(0, window):
        value = float("{0:.3f}".format(abs(sign_cplx[h])))
        sign_cplx_abs.append(value)
    print("")
    return sign_cplx, sign_cplx_abs


def mult_signal_on_constant(signal=None,
                            window=None,
                            harmonic=None,
                            samp_rate=None,
                            print_events=True):
    if print_events:
        print("\n[def] Shift in the signal " + str(harmonic) +
              " harmonica to the center, in the band " + str(samp_rate) + " Hz")
    start_time = time.time()

    signal_shift = []
    carrier = ((window - harmonic) * samp_rate) / window

    step_t = float(1 / samp_rate)
    for i in range(0, window):
        cos_ = float(math.cos((2 * np.pi * carrier * step_t * i)))
        sin_ = float(math.sin((2 * np.pi * carrier * step_t * i)))
        point = complex(cos_, sin_)
        new_sample = point * signal[i]
        signal_shift.append(new_sample)

    signal_shift_abs = []
    for sps in signal_shift:
        z = complex((sps.real), (sps.imag))
        value = float("{0:.3f}".format(abs(z)))
        signal_shift_abs.append(value)
    if print_events:
        print("--- %s seconds ---" % (time.time() - start_time))
    return signal_shift, signal_shift_abs


def make_filtration(signal=None,
                    window=None,
                    harmonic=None,
                    samp_rate=None,
                    filtration_width=None,
                    on_center=None,
                    print_event=True):
    if print_event:
        print("\n [def] Signal filtering")
        print(" Filter width: " + str(filtration_width) + ", central harmonica: " +
              str(harmonic) + ", полоса: " + str(samp_rate) + " Hz")
    start_time = time.time()

    signal_before_filtration, d = mult_signal_on_constant(signal=signal, window=window, harmonic=harmonic,
                                                          samp_rate=samp_rate, print_events=False)
    coeff_f_spectrum_cplx, coeff_f_spectrum_abs = fast_ft(signal=signal_before_filtration, window=window,
                                                          output=None, print_events=False)

    for i in range(int(filtration_width / 2), int(window - filtration_width / 2)):
        coeff_f_spectrum_cplx[i] = complex(0)

    signal_cplx_filtered, signal_cplx_abs_filtered = fast_ft_reverse(coefficients=coeff_f_spectrum_cplx,
                                                                     window=window, print_events=False)

    if on_center:
        signal_after_centred, signal_after_centred_abs = mult_signal_on_constant(signal=signal_cplx_filtered,
                                                                                 window=window,
                                                                                 harmonic=int(window / 2),
                                                                                 samp_rate=samp_rate,
                                                                                 print_events=False)
        if print_event:
            print("--- %s seconds ---" % (time.time() - start_time))
        return signal_after_centred, signal_after_centred_abs
    else:
        if print_event:
            print("--- %s seconds ---" % (time.time() - start_time))
        return signal_cplx_filtered, signal_cplx_abs_filtered


def get_shifted_window(signal=None,
                       shift=None,
                       window=None,
                       half_shift=False):
    sign_cplx = []
    if half_shift:
        if shift == 0:
            for i in range(0, window):
                sign_cplx.append(complex(signal[i * 2], signal[i * 2 + 1]))
        else:
            for i in range(int(shift * (window / 2)), int((shift * (window / 2)) + window)):
                sign_cplx.append(complex(signal[i * 2], signal[i * 2 + 1]))
    else:
        for i in range(int(shift * window), int((shift + 1) * window)):
            sign_cplx.append(complex(signal[i * 2], signal[i * 2 + 1]))

    return sign_cplx


def save_to_file_filtered_signal(signal=None,
                                 window=None,
                                 samp_rate=None,
                                 filtr_width=None,
                                 on_center=False,
                                 harmonic=None,
                                 count_windows=None,
                                 output=None,
                                 print_event=True,
                                 safe_format="float32"):
    if print_event:
        print("\n[def] Save the new file to " + output + "_" + safe_format + ".bin" + ", windows: " + str(count_windows))

    start_time = time.time()
    filtered_signal_parted = []

    for i in range(0, count_windows):
        part_of_signal = get_shifted_window(signal=signal,
                                            window=window,
                                            shift=i)

        filtered_part_signal, filtered_part_signal_abs = make_filtration(signal=part_of_signal,
                                                                         window=window,
                                                                         harmonic=harmonic,
                                                                         samp_rate=samp_rate,
                                                                         filtration_width=filtr_width,
                                                                         on_center=on_center,
                                                                         print_event=False)

        filtered_signal_parted.append(filtered_part_signal)

    filtered_signal_few_window = []
    for j in filtered_signal_parted:
        for k in j:
            filtered_signal_few_window.append(k.real)
            filtered_signal_few_window.append(k.imag)

    with open(output + "_" + safe_format + ".bin", "wb") as pcm_file:
        for samples in filtered_signal_few_window:
            if safe_format == "float32":
                bytearray_sample = np.float32(samples / window)
            elif safe_format == "float16":
                bytearray_sample = np.float16(samples / window)
            else:
                bytearray_sample = bytearray(struct.pack("f", samples / window))
            pcm_file.write(bytearray_sample)
    if print_event:
        print(" Recorded", len(filtered_signal_few_window), "samples. Windows:", count_windows)
        print("--- %s seconds ---" % (time.time() - start_time))


def save_to_file_filtered_signal_window(signal=None,
                                        window=None,
                                        samp_rate=None,
                                        filtr_width=None,
                                        on_center=False,
                                        harmonic=None,
                                        count_windows=None,
                                        output=None,
                                        wind_var=None,
                                        print_event=True,
                                        safe_format="float32"):
    start_time = time.time()

    window_string = ""
    if wind_var == 1:
        window_string = "sin"
    elif wind_var == 2:
        window_string = "hann"
    elif wind_var == 3:
        window_string = "hemming"
    elif wind_var == 4:
        window_string = "blkmn-harris"
    elif wind_var == 5:
        window_string = "blkmn-natoll"

    function_windowed_array = []
    for i in range(0, window):
        if wind_var == 1:
            function_windowed_array.append(math.sin((np.pi * i) / (window - 1)))
        elif wind_var == 2:
            function_windowed_array.append(0.5 * (1 - math.cos((2 * np.pi * i) / (window - 1))))
        elif wind_var == 3:
            function_windowed_array.append(0.53836 - (0.46164 * math.cos((2 * np.pi * i) / (window - 1))))
        elif wind_var == 4:
            a0 = 0.35875
            a1 = 0.48829
            a2 = 0.14128
            a3 = 0.01168
            function_windowed_array.append(a0 - a1 * math.cos((2 * np.pi * i) / (window - 1))
                                              + a2 * math.cos((4 * np.pi * i) / (window - 1))
                                              - a3 * math.cos((6 * np.pi * i) / (window - 1)))
        elif wind_var == 5:
            a0 = 0.3635819
            a1 = 0.4891775
            a2 = 0.1365995
            a3 = 0.0106411
            function_windowed_array.append(a0 - a1 * math.cos((2 * np.pi * i) / (window - 1))
                                              + a2 * math.cos((4 * np.pi * i) / (window - 1))
                                              - a3 * math.cos((6 * np.pi * i) / (window - 1)))

    if print_event:
        print("\n[def] Save the new file to " + output + "_" + window_string + "_" + safe_format + ".bin")

    fst_of_signal = get_shifted_window(signal=signal,
                                       window=window,
                                       shift=0,
                                       half_shift=False)

    filtered_part_signal, filtered_part_signal_abs = make_filtration(signal=fst_of_signal,
                                                                     window=window,
                                                                     harmonic=harmonic,
                                                                     samp_rate=samp_rate,
                                                                     filtration_width=filtr_width,
                                                                     on_center=on_center,
                                                                     print_event=False)

    fst_multilpied_array_window_with_filtered = []
    for j in range(0, window):
        current_value = function_windowed_array[j] * filtered_part_signal[j]
        fst_multilpied_array_window_with_filtered.append(current_value)

    fst_half_array = fst_multilpied_array_window_with_filtered[:int(window / 2)]

    with open(output + "_" + window_string + "_" + safe_format + ".bin", "wb") as pcm_file:
        for samples in fst_half_array:
            real_sample = (samples.real) / window
            imag_sample = (samples.imag) / window
            if safe_format == "float32":
                bytearray_real_sample = np.float32(real_sample)
                bytearray_imag_sample = np.float32(imag_sample)
            elif safe_format == "float16":
                bytearray_real_sample = np.float16(real_sample)
                bytearray_imag_sample = np.float16(imag_sample)
            else:
                bytearray_real_sample = bytearray(struct.pack("f", real_sample))
                bytearray_imag_sample = bytearray(struct.pack("f", imag_sample))
            pcm_file.write(bytearray_real_sample)
            pcm_file.write(bytearray_imag_sample)

    snd_half_array = fst_multilpied_array_window_with_filtered[int(window / 2):]

    for i in range(1, count_windows):
        part_of_signal = get_shifted_window(signal=signal,
                                            window=window,
                                            shift=i,
                                            half_shift=True)

        filtered_part_signal, filtered_part_signal_abs = make_filtration(signal=part_of_signal,
                                                                         window=window,
                                                                         harmonic=harmonic,
                                                                         samp_rate=samp_rate,
                                                                         filtration_width=filtr_width,
                                                                         on_center=on_center,
                                                                         print_event=False)

        multilpied_array_windowed_with_filtered = []
        for j in range(0, window):
            current_value = function_windowed_array[j] * filtered_part_signal[j]
            multilpied_array_windowed_with_filtered.append(current_value)

        multilpied_array_windowed_with_filtered_inverse = multilpied_array_windowed_with_filtered[::-1]
        half_array_first_part = multilpied_array_windowed_with_filtered_inverse[:int(window / 2)]

        new_array_to_file = []
        for i in range(0, int(window / 2)):
            new_array_to_file.append(0)

        for i in range(0, int(window / 2)):
            new_array_to_file[i] = (half_array_first_part[i] + snd_half_array[i])


        with open(output + "_" + window_string + "_" + safe_format + ".bin", "ab") as pcm_file:
            for samples in new_array_to_file:
                real_sample = samples.real / window
                imag_sample = samples.imag / window
                if safe_format == "float32":
                    bytearray_real_sample = np.float32(real_sample)
                    bytearray_imag_sample = np.float32(imag_sample)
                elif safe_format == "float16":
                    bytearray_real_sample = np.float16(real_sample)
                    bytearray_imag_sample = np.float16(imag_sample)
                else:
                    bytearray_real_sample = bytearray(struct.pack("f", real_sample))
                    bytearray_imag_sample = bytearray(struct.pack("f", imag_sample))
                pcm_file.write(bytearray_real_sample)
                pcm_file.write(bytearray_imag_sample)

        snd_half_array = multilpied_array_windowed_with_filtered_inverse[int(window / 2):]


    if print_event:
        print("--- %s seconds ---" % (time.time() - start_time))
