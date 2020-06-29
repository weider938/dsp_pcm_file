# Ctrl + Alt + L --> PEP8
# author: weider9385@gmail.com

from dsp_lib import *
import sys

if __name__ == "__main__":
    file_input = None
    my_signal = None
    chunks_cnt = None
    deep_FT = None
    mode = None

    try:
        file_input = sys.argv[1]
        my_signal, chunks_cnt = open_bin_file(f_name=file_input)
        print(" File input: " + file_input)
    except Exception as e:
        try:
            my_signal, chunks_cnt = open_bin_file(f_name="333.pcm")
            print(" File input: 333.pcm")
        except Exception as e:
            try:
                print(" Error: the source file was not written as an argument after the executable.")
                print(" Examples: \"python interactive.py 222.pcm\" или DSP_Console.exe 444.pcm")
                print(" Error reading file: ", str(e))
                print(" Enter a file name:", end="")
                new_file = input()
                my_signal, chunks_cnt = open_bin_file(f_name=new_file)
                print(" File input: " + new_file)
            except:
                print("")
                print(" Shutdown.")
                print("")
                sys.exit()
    try:
        deep_FT = int(input("\n Enter window size for FT (for FFT: 128,..., 4096, 8192..., 262144 (2^18)): "))
    except:
        try:
            deep_FT = int(input("\n Enter window size for FT (for FFT: 128,..., 4096, 8192..., 262144 (2^18)): "))
        except:
            deep_FT = 8192
            print(" The offer to enter the PF window was ignored, the value 8192 was selected")
    
    print("\n Select an operating mode:")
    print("  0 - Write spectrum coeffs to .txt with Discrete FT;")
    print("  1 - Write spectrum coeffs to .txt with Fast FT;")
    print("  2 - Write spectrum coeffs to .txt with a carrier frequency offset of 0;")
    print("  3 - Write spectrum coeffs of the filtered signal (in the center);")
    print("  4 - Signal Filtering (no windows) and save to binary;")
    print("  5 - Signal Filtering using window functions and save to binary;")
    print("  6 - Write coeffs to .txt squared spectrum;")
    print("  7 - Write coeffs to .txt of 4-degree spectrum construction")
    print("  8 - Show MENU;")
    print("  9 - Exit.")

    while True:

        try:
            mode = int(input("\n Enter the digits [0..9] (menu ->8): "))
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 0:
                coeff_spectrum_abs_slow = slow_ft(signal=my_signal, window=deep_FT, output="ft_slow_spectrum.txt")
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 1:
                signal_cplx, signal_cplx_abs = return_samples(signal=my_signal,
                                                              window=deep_FT)
                coeff_spectrum_cplx, coeff_spectrum_abs = fast_ft(signal=signal_cplx,
                                                                  window=deep_FT,
                                                                  output="ft_fast_spectrum.txt")
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 2:
                harmonic = int(input(" Enter the number of the harmonic to be shifted (e.g. 6744): "))
                signal_cplx, signal_cplx_abs = return_samples(signal=my_signal,
                                                              window=deep_FT)
                sign_with_shift, sign_with_shift_abs = mult_signal_on_constant(signal=signal_cplx,
                                                                               window=deep_FT,
                                                                               harmonic=harmonic,
                                                                               samp_rate=2 * (10 ** 8))
                coeff_spectrum_cplx_shift, coeff_spectrum_abs_shift = fast_ft(signal=sign_with_shift,
                                                                              window=deep_FT,
                                                                              output="shifted_signal_spectrum.txt",
                                                                              print_events=False)
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 3:
                signal_cplx, signal_cplx_abs = return_samples(signal=my_signal,
                                                              window=deep_FT)
                harmonic = int(input(" Enter the number of the harmonic to be shifted (e.g. 6744): "))
                width = int(input(" Enter the filter width (in harmonics, e.g. 1000): "))
                filtered_signal, filtered_signal_abs = make_filtration(signal=signal_cplx,
                                                                       window=deep_FT,
                                                                       harmonic=harmonic,
                                                                       on_center=True,
                                                                       samp_rate=2 * (10 ** 8),
                                                                       filtration_width=width)

                coeff_spectrum_cplx_filtered, coeff_spectrum_abs_filtered = fast_ft(signal=filtered_signal,
                                                                                    window=deep_FT,
                                                                                    output="filtered_signal_spectrum_on_center.txt",
                                                                                    print_events=False)
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 4:
                harmonic = int(input(" Enter the number of the harmonic to be shifted (e.g. 6744): "))
                width = int(input(" Enter the filter width (in harmonics, e.g. 1000): "))
                count_win = int(input(" Enter the number of windows:"))

                save_to_file_filtered_signal(signal=my_signal,
                                             window=deep_FT,
                                             samp_rate=2 * (10 ** 8),
                                             filtr_width=width,
                                             on_center=False,
                                             harmonic=harmonic,
                                             count_windows=count_win,
                                             output="filtered_without_window",
                                             safe_format="float32")  # float32 || float
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 5:
                window_variable = int(input(" Выберите окно \n  1 - sinus (показывает лучшие свойства),"
                                            " \n  2 - Hann (Hanning),"
                                            " \n  3 - Hemming,"
                                            " \n  4 - Blackman–Harris,"
                                            " \n  5 - Blackman–Natoll: "))

                harmonic = int(input(" Enter the number of the harmonic to be shifted (e.g. 6744): "))
                width = int(input(" Enter the filter width (in harmonics, 1000 for example): "))
                count_win = int(input(" Enter the number of windows:"))

                save_to_file_filtered_signal_window(signal=my_signal,
                                                     window=deep_FT,
                                                     samp_rate=2 * (10 ** 8),
                                                     filtr_width=width,
                                                     on_center=False,
                                                     harmonic=harmonic,
                                                     count_windows=count_win*2,
                                                     output="filtered_with_window",
                                                     wind_var = window_variable,
                                                     safe_format="float32")  # float16 || float32 || float
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 6:
                harmonic = int(input(" Enter the number of the harmonic to be shifted (e.g. 6744): "))
                width = int(input(" Enter the filter width (in harmonics, e.g. 1000): "))
                signal_cplx, signal_cplx_abs = return_samples(signal=my_signal,
                                                              window=deep_FT)

                filtered_signal, filtered_signal_abs = make_filtration(signal=signal_cplx,
                                                                       window=deep_FT,
                                                                       harmonic=harmonic,
                                                                       on_center=False,
                                                                       samp_rate=2 * (10 ** 8),
                                                                       filtration_width=width)


                signal_after_centred, signal_after_centred_abs = mult_signal_on_constant(signal=filtered_signal,
                                                                                         window=deep_FT,
                                                                                         harmonic=int(deep_FT / 2),
                                                                                         samp_rate=2 * (10 ** 8),
                                                                                         print_events=False)

                quadr_signal = []
                for i in range(0, deep_FT):
                    quadr_signal.append(
                        signal_after_centred[i] * signal_after_centred[i])

                signal_after_centred5, signal_after_centred_abs5 = mult_signal_on_constant(signal=quadr_signal,
                                                                                           window=deep_FT,
                                                                                           harmonic=int(deep_FT / 2),
                                                                                           samp_rate=2 * (10 ** 8),
                                                                                           print_events=False)

                coeff_spectrum_cplx6, coeff_spectrum_abs6 = fast_ft(signal=signal_after_centred5,
                                                                    window=deep_FT,
                                                                    output="spectrum_2_step.txt")
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 7:
                harmonic = int(input(" Enter the number of the harmonic to be shifted (e.g. 6744): "))
                width = int(input(" Enter the filter width (in harmonics, e.g. 1000): "))
                signal_cplx, signal_cplx_abs = return_samples(signal=my_signal,
                                                              window=deep_FT)

                filtered_signal, filtered_signal_abs = make_filtration(signal=signal_cplx,
                                                                       window=deep_FT,
                                                                       harmonic=harmonic,
                                                                       on_center=False,
                                                                       samp_rate=2 * (10 ** 8),
                                                                       filtration_width=width)


                signal_after_centred, signal_after_centred_abs = mult_signal_on_constant(signal=filtered_signal,
                                                                                         window=deep_FT,
                                                                                         harmonic=int(deep_FT / 2),
                                                                                         samp_rate=2 * (10 ** 8),
                                                                                         print_events=False)

                quadr_signal = []
                for i in range(0, deep_FT):
                    quadr_signal.append(signal_after_centred[i]
                                        * signal_after_centred[i]
                                        * signal_after_centred[i]
                                        * signal_after_centred[i])


                signal_after_centred5, signal_after_centred_abs5 = mult_signal_on_constant(signal=quadr_signal,
                                                                                           window=deep_FT,
                                                                                           harmonic=int(deep_FT / 2),
                                                                                           samp_rate=2 * (10 ** 8),
                                                                                           print_events=False)

                coeff_spectrum_cplx6, coeff_spectrum_abs6 = fast_ft(signal=signal_after_centred5,
                                                                    window=deep_FT,
                                                                    output="spectrum_4_step.txt")
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 8:
    print("\n Select an operating mode:")
    print("  0 - Write spectrum coeffs to .txt with Discrete FT;")
    print("  1 - Write spectrum coeffs to .txt with Fast FT;")
    print("  2 - Write spectrum coeffs to .txt with a carrier frequency offset of 0;")
    print("  3 - Write spectrum coeffs of the filtered signal (in the center);")
    print("  4 - Signal Filtering (no windows) and save to binary;")
    print("  5 - Signal Filtering using window functions and save to binary;")
    print("  6 - Write coeffs to .txt squared spectrum;")
    print("  7 - Write coeffs to .txt of 4-degree spectrum construction")
    print("  8 - Show MENU;")
    print("  9 - Exit.")
        except Exception as e:
            print(" Error: " + str(e))

        try:
            if mode == 9:
                print("")
                print(" Shutdown.")
                print("")
                sys.exit()
        except Exception as e:
            print(" Error: " + str(e))