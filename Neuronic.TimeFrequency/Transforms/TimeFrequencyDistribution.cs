﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using Neuronic.TimeFrequency.Kernels;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Time-Frequency distribution (TFD).
    /// </summary>
    public class TimeFrequencyDistribution : IBilinearTimeFrequencyRepresentation
    {
        private readonly double[,] _values;
        private readonly double[] _frequencies;

        internal TimeFrequencyDistribution(double[,] values, double[] frequencies, 
            double startTime, double samplingPeriod)
        {
            _values = values;
            _frequencies = frequencies;
            SamplingPeriod = samplingPeriod;
            StartTime = startTime;
        }

        private static Complex[] GetAnalytic(IReadOnlySignal<double> signal)
        {
            var z = new Signal<Complex>(new Complex[2 * signal.Count], signal.Start, signal.SamplingRate);
            for (int i = 0; i < signal.Count; i++)
                z[i] = signal[i];
            
            z.HilbertTransform();

            for (int i = 0; i < signal.Count; i++)
                z[signal.Count + i] = 0;

            return z.Samples;
        }

        private static List<Complex[]> GetTL(int n, int n2, Complex[] z)
        {
            var kTrans = new List<Complex[]>(2 * n);
            for (int j = 0; j < (n + 1) / 2; j++)
            {
                var column = new Complex[n];
                for (int i = 0; i < n; i++)
                {
                    var i1 = (i + j) % n2;
                    var i2 = (i - j) % n2;
                    if (i2 < 0) i2 += n2;

                    column[i] = z[i1] * Complex.Conjugate(z[i2]);
                }

                kTrans.Add(column);
                column = new Complex[n];
                for (int i = 0; i < n; i++)
                {
                    var i2 = (i - j) % n2;
                    if (i2 < 0) i2 += n2;
                    var i3 = (i + j + 1) % n2;

                    column[i] = z[i3] * Complex.Conjugate(z[i2]);
                }

                kTrans.Add(column);
            }

            return kTrans;
        }

        private static void Convolve(int n, int n2, List<Complex[]> kTrans, double[,] kernel, ParallelOptions options)
        {
            while (kTrans.Count < 2 * n)
                kTrans.Add(new Complex[n]);

            Parallel.For(0, n + 1, options, j =>
            {
                var column = kTrans[j];
                // fft
                column.FFT();
                // conv
                for (int i = 0; i < n; i++)
                    column[i] *= kernel[i, j];
                // ifft
                column.IFFT();
            });

            for (int i = 0; i < n; i++)
            {
                for (int j = 1; j < (n + 2) / 2; j++)
                    kTrans[n2 - 2 * j][i] = Complex.Conjugate(kTrans[2 * j][i]);
                for (int j = 0; j < (n + 2) / 2; j++)
                    kTrans[n2 - 2 * j - 1][i] = Complex.Conjugate(kTrans[2 * j + 1][i]);
            }
        }

        private static double[,] ToTFDomain(int freq, int time, List<Complex[]> kTrans, ParallelOptions options)
        {
            var tfd = new double[time, freq];

            Parallel.For(0, freq, options,
                () => new {buffer = new Complex[freq]},
                (i, _, state) =>
                {
                    var buffer = state.buffer;

                    // Even rows
                    for (int j = 0; j < freq; j++)
                        buffer[j] = kTrans[2 * j][i];
                    buffer.FFT();
                    for (int j = 0; j < freq; j++)
                        tfd[2 * i, j] = (buffer[j] / time).Real;

                    for (int j = 0; j < freq; j++)
                        buffer[j] = kTrans[2 * j + 1][i];

                    // Odd rows
                    tfd[2 * i + 1, 0] = (Enumerable.Range(0, freq).Aggregate(Complex.Zero, (total, j) => total + buffer[j]) / time).Real;

                    var nFreq = (freq + 1) / 2;
                    buffer[0] = buffer[0].Imaginary;
                    for (int k = 1; k <= nFreq; k++)
                        buffer[k] = (buffer[k] - Complex.Conjugate(buffer[freq - k])) / new Complex(0, 2);
                    for (int ke = nFreq + 1; ke < freq; ke++)
                        buffer[ke] = Complex.Conjugate(buffer[freq - ke]);


                    buffer.FFT();

                    for (int j = 1; j < freq; j++)
                    {
                        var a = Math.Cos((Math.PI / freq) * j);
                        var b = Math.Sin((Math.PI / freq) * j);
                        var c = (a * a + b * b) / b;

                        tfd[2 * i + 1, j] = (buffer[j] * c / time).Real;
                    }

                    return state;
                },
                _ => { });

            return tfd;
        }

        /// <summary>
        /// Estimates the TFD of the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="kernel">The kernel.</param>
        /// <param name="options">The options for parallelization.</param>
        /// <returns>The estimated TFD.</returns>
        /// <exception cref="ArgumentNullException">Thrown when either the signal or the kernel are <c>null</c>.</exception>
        /// <remarks>
        /// This algorithm is based on the one proposed by J.M. O' Toole and B. Boashash on the article
        /// "Fast and memory-efficient algorithms for computing quadratic time–frequency distributions",
        /// Applied and Computational Harmonic Analysis, vol. 35, no. 2, pp. 350–358, 2013.
        /// </remarks>
        public static TimeFrequencyDistribution Estimate(IReadOnlySignal<double> signal, DopplerLagKernel kernel, ParallelOptions options = null)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            if (kernel == null) throw new ArgumentNullException(nameof(kernel));
            options = options ?? new ParallelOptions();

            var z = GetAnalytic(signal);
            var n2 = z.Length;
            var n = n2 / 2;

            var g = kernel.Evaluate(n);
            
            var kTrans = GetTL(n, n2, z);

            Convolve(n, n2, kTrans, g, options);

            var tfd = ToTFDomain(n, n2, kTrans, options);

            var frequencies = new double[tfd.GetLength(1)];
            for (int i = 0; i < frequencies.Length; i++)
                frequencies[i] = (i * 0.5 * signal.SamplingRate) / (frequencies.Length - 1);
            return new TimeFrequencyDistribution(tfd, frequencies, signal.Start, signal.SamplingPeriod / 2);
        }

        /// <summary>
        /// Estimates the TFD of the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="kernel">The kernel.</param>
        /// <param name="options">The options for parallelization.</param>
        /// <returns>The estimated TFD.</returns>
        /// <exception cref="ArgumentNullException">Thrown when either the signal or the kernel are <c>null</c>.</exception>
        /// <remarks>
        /// This algorithm is based on the one proposed by J.M. O' Toole and B. Boashash on the article
        /// "Fast and memory-efficient algorithms for computing quadratic time–frequency distributions",
        /// Applied and Computational Harmonic Analysis, vol. 35, no. 2, pp. 350–358, 2013.
        /// </remarks>
        public static TimeFrequencyDistribution Estimate(IReadOnlySignal<float> signal, DopplerLagKernel kernel, ParallelOptions options = null)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            return Estimate(signal.Map(x => (double) x), kernel, options);
        }

        /// <summary>
        /// Gets the offset of the first sample in the time domain.
        /// </summary>
        public double StartTime { get; }

        /// <summary>
        /// Gets the sampling period.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the amount of samples.
        /// </summary>
        public int SampleCount => _values.GetLength(0);

        /// <summary>
        /// Gets the amount of samples in the frequency domain.
        /// </summary>
        public int FrequencyCount => _values.GetLength(1);

        /// <summary>
        /// Gets the computed frequencies.
        /// </summary>
        public IEnumerable<double> Frequencies => _frequencies.AsReadOnly();

        double IBilinearTimeFrequencyRepresentation.this[int offset, int frequencyIndex] =>
            this[offset, frequencyIndex];

        /// <summary>
        /// Gets the TFD value for the specified frequency and offset.
        /// </summary>
        /// <param name="offset">The offset.</param>
        /// <param name="freqIndex">The frequency index.</param>
        /// <returns>The TFD value.</returns>
        public double this[int offset, int freqIndex] => _values[offset, freqIndex];

        double ITimeFrequencyRepresentation.this[int offset, double frequency]
        {
            get
            {
                var freqIndex = Array.BinarySearch(_frequencies, frequency);
                if (freqIndex < 0)
                    freqIndex = ~freqIndex;
                return this[offset, freqIndex];
            }
        }

        /// <summary>
        /// Gets the TFD value for the specified frequency and offset.
        /// </summary>
        /// <param name="delay">The delay.</param>
        /// <param name="frequency">The frequency.</param>
        /// <returns>The TFD value.</returns>
        public double this[double delay, double frequency]
        {
            get
            {
                var offset = this.FindClosestOffsetOfTime(delay);
                var freqIndex = this.FindClosestIndexOfFrequency(frequency);
                return this[offset, freqIndex];
            }
        }
    }
}