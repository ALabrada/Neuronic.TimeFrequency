using System;
using System.Collections.Generic;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;

namespace Neuronic.TimeFrequency
{
    class ShortTimeFourierTransform
    {
        private Complex[,] _values;
        private double[] _frequencies;

        public ShortTimeFourierTransform(float[] samples, TimeSpan samplingPeriod, double[] window)
        {
            SamplingPeriod = samplingPeriod;
            var wFreq = 2 * Math.PI / samplingPeriod.TotalSeconds;
            var count = samples.Length / 2;
            _frequencies = new double[count + 1];
            for (int i = 0; i <= count; i++)
                _frequencies[i] = wFreq * i / samples.Length;

            _values = new Complex[samples.Length, count + 1];
            var windowed = new Complex[samples.Length];

            for (int offset = 0; offset < samples.Length; offset++)
            {
                for (int i = 0; i < window.Length; i++)
                {
                    var index = offset + i - window.Length / 2;
                    if (index >= 0 && index < samples.Length)
                        windowed[index] = window[i] * samples[index];
                }

                FourierTransform2.DFT(windowed, FourierTransform.Direction.Forward);

                for (int i = 0; i <= count; i++)
                    _values[offset, i] = samplingPeriod.TotalSeconds * windowed[i];

                Array.Clear(windowed, 0, windowed.Length);
            }
        }

        public ShortTimeFourierTransform(float[] samples, TimeSpan samplingPeriod, int windowRadius)
            : this(samples, samplingPeriod, GaussWindow(windowRadius, samplingPeriod))
        {
        }

        public static double[] GaussWindow(int radius, TimeSpan samplingPeriod)
        {
            var window = new double[2 * radius + 1];
            var factor = 1 / Math.Pow(Math.PI, 0.25);
            for (int i = 0; i <= radius; i++)
            {
                var t = i * samplingPeriod.TotalSeconds;
                window[i] = window[2 * radius - i] = factor * Math.Exp(-t * t / 2);
            }
            return window;
        }

        public TimeSpan SamplingPeriod { get; }

        public IEnumerable<double> Frequencies => _frequencies;

        public Complex this[int offset, int freqIndex] => _values[offset, freqIndex];

        public Complex this[double delay, double frequency]
        {
            get
            {
                var offset = (int)Math.Round(delay / SamplingPeriod.TotalSeconds);
                var freqIndex = Array.BinarySearch(_frequencies, frequency);
                if (freqIndex < 0)
                    freqIndex = ~freqIndex;
                return this[offset, freqIndex];
            }
        }
    }
}
