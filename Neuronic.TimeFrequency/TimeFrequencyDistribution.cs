using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;
using Neuronic.TimeFrequency.Kernels;

namespace Neuronic.TimeFrequency
{
    public class TimeFrequencyDistribution
    {
        private readonly double[,] _values;
        private readonly double[] _frequencies;

        private TimeFrequencyDistribution(double[,] values, double[] frequencies, double samplingPeriod)
        {
            _values = values;
            _frequencies = frequencies;
            SamplingPeriod = samplingPeriod;
        }

        private static Complex[] GetAnalytic(IReadOnlySignal<double> signal)
        {
            var z = new Complex[2 * signal.Count];
            for (int i = 0; i < signal.Count; i++)
                z[i] = signal[i];
            FourierTransform2.FFT(z, FourierTransform.Direction.Forward);

            for (int i = 1; i < signal.Count; i++)
            {
                z[i] *= 2d;
                z[signal.Count + i] = 0;
            }

            FourierTransform2.FFT(z, FourierTransform.Direction.Backward);
            for (int i = 0; i < signal.Count; i++)
                z[signal.Count + i] = 0;

            return z;
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

        private static void Convolve(int n, int n2, List<Complex[]> kTrans, double[,] kernel)
        {
            while (kTrans.Count < 2 * n)
                kTrans.Add(new Complex[n]);

            for (int j = 0; j < n + 1; j++)
            {
                var column = kTrans[j];
                // fft
                FourierTransform2.FFT(column, FourierTransform.Direction.Forward);
                // conv
                for (int i = 0; i < n; i++)
                    column[i] *= kernel[i, j];
                // ifft
                FourierTransform2.FFT(column, FourierTransform.Direction.Backward);
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 1; j < (n + 2) / 2; j++)
                    kTrans[n2 - 2 * j][i] = Complex.Conjugate(kTrans[2 * j][i]);
                for (int j = 0; j < (n + 2) / 2; j++)
                    kTrans[n2 - 2 * j - 1][i] = Complex.Conjugate(kTrans[2 * j + 1][i]);
            }
        }

        private static double[,] ToTFDomain(int freq, int time, List<Complex[]> kTrans)
        {
            var tfd = new double[time, freq];
            var buffer = new Complex[freq];
            for (int i = 0; i < freq; i++)
            {
                for (int j = 0; j < freq; j++)
                    buffer[j] = kTrans[2 * j][i];
                FourierTransform2.FFT(buffer, FourierTransform.Direction.Forward);
                for (int j = 0; j < freq; j++)
                    tfd[2 * i, j] = (buffer[j] / time).Real;
            }

            for (int i = 0; i < freq; i++)
            {
                for (int j = 0; j < freq; j++)
                    buffer[j] = kTrans[2 * j + 1][i];

                tfd[2 * i + 1, 0] = (Enumerable.Range(0, freq).Aggregate(Complex.Zero, (total, j) => total + buffer[j]) / time).Real;

                var nFreq = (freq + 1) / 2;
                buffer[0] = buffer[0].Imaginary;
                for (int k = 1; k <= nFreq; k++)
                    buffer[k] = (buffer[k] - Complex.Conjugate(buffer[freq - k])) / new Complex(0, 2);
                for (int ke = nFreq + 1; ke < freq; ke++)
                    buffer[ke] = Complex.Conjugate(buffer[freq - ke]);
                

                FourierTransform2.FFT(buffer, FourierTransform.Direction.Forward);

                for (int j = 1; j < freq; j++)
                {
                    var a = Math.Cos((Math.PI / freq) * j);
                    var b = Math.Sin((Math.PI / freq) * j);
                    var c = (a * a + b * b) / b;

                    tfd[2 * i + 1, j] = (buffer[j] * c / time).Real;
                }
            }

            return tfd;
        }

        public static TimeFrequencyDistribution Estimate(IReadOnlySignal<double> signal, NonSeparableKernel kernel)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            if (kernel == null) throw new ArgumentNullException(nameof(kernel));

            var z = GetAnalytic(signal);
            var n2 = z.Length;
            var n = n2 / 2;

            var g = kernel.Evaluate(n);
            
            var kTrans = GetTL(n, n2, z);

            Convolve(n, n2, kTrans, g);

            var tfd = ToTFDomain(n, n2, kTrans);

            var frequencies = new double[tfd.GetLength(1)];
            for (int i = 0; i < frequencies.Length; i++)
                frequencies[i] = (i * 0.5 * signal.SamplingRate) / (frequencies.Length - 1);
            return new TimeFrequencyDistribution(tfd, frequencies, signal.SamplingPeriod / 2);
        }

        public static TimeFrequencyDistribution Estimate(IReadOnlySignal<float> signal, NonSeparableKernel kernel)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            return Estimate(signal.Map(x => (double) x), kernel);
        }

        public double SamplingPeriod { get; }

        public int Samples => _values.GetLength(0);

        public ReadOnlyCollection<double> Frequencies => Array.AsReadOnly(_frequencies);

        public Complex this[int offset, int freqIndex] => _values[offset, freqIndex];

        public Complex this[double delay, double frequency]
        {
            get
            {
                var offset = (int)Math.Round(delay / SamplingPeriod);
                var freqIndex = Array.BinarySearch(_frequencies, frequency);
                if (freqIndex < 0)
                    freqIndex = ~freqIndex;
                return this[offset, freqIndex];
            }
        }

        public IEnumerable<double> EnumerateFrequency(int index) =>
            Enumerable.Range(0, _values.GetLength(0)).Select(j => _values[j, index]);

        public IEnumerable<double> EnumerateOffset(int index) =>
            Enumerable.Range(0, _values.GetLength(1)).Select(i => _values[index, i]);
    }
}