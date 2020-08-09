using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;

namespace Neuronic.TimeFrequency
{
    class ContinuousWaveletTransform: IEnumerable<Complex>
    {
        private Complex[,] _values;
        private double[] _scales;

        private ContinuousWaveletTransform(Complex[,] values, TimeSpan samplingPeriod, IWavelet wavelet, IEnumerable<double> scales)
        {
            _values = values;
            SamplingPeriod = samplingPeriod;
            Wavelet = wavelet;
            _scales = scales.ToArray();
        }

        public static ContinuousWaveletTransform EstimateUsingFFT(float[] samples, TimeSpan samplingPeriod, IWavelet wavelet, IEnumerable<double> scales)
        {
            if (samples is null)            
                throw new ArgumentNullException(nameof(samples));  
            
            scales = scales ?? Enumerable.Range(1, samples.Length).Select(i => samplingPeriod.TotalSeconds * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Haar;

            var general = new Complex[samples.Length];
            for (int offset = 0; offset < samples.Length; offset++)
                general[offset] = samples[offset];
            FourierTransform2.DFT(general, FourierTransform.Direction.Forward);


            var values = new Complex[samples.Length, scaleArray.Length];
            var scaled = new Complex[samples.Length];

            for (int scale = 0; scale < scaleArray.Length; scale++)
            {
                var dt = (samples.Length - 1) * samplingPeriod.TotalSeconds / scaleArray[scale];
                wavelet.Evaluate(-dt, 0, scaled, 0, scaled.Length);
                Array.Reverse(scaled, 0, scaled.Length);

                for (int offset = 0; offset < samples.Length; offset++)
                    scaled[offset] = Complex.Conjugate(scaled[offset]);

                FourierTransform2.DFT(scaled, FourierTransform.Direction.Forward);

                var factor = 1d / Complex.Sqrt(wavelet.Energy * Math.Abs(scaleArray[scale]));
                for (int offset = 0; offset < samples.Length; offset++)
                    scaled[offset] = factor * scaled[offset] * general[offset];

                FourierTransform2.DFT(scaled, FourierTransform.Direction.Backward);

                for (int offset = 0; offset < samples.Length; offset++)
                    values[offset, scale] = scaled[offset];
            }

            return new ContinuousWaveletTransform(values, samplingPeriod, wavelet, scaleArray);
        }

        public TimeSpan SamplingPeriod { get; }

        public IWavelet Wavelet { get; }

        public IReadOnlyList<double> Scales => _scales;

        public Complex this[int offset, int scaleIndex] => _values[offset, scaleIndex];

        public Complex this[double delay, double scale]
        {
            get
            {
                var offset = (int)Math.Round(delay / SamplingPeriod.TotalSeconds);
                var scaleIndex = Array.BinarySearch(_scales, scale);
                if (scaleIndex < 0)
                    scaleIndex = ~scaleIndex;
                return this[offset, scaleIndex];
            }
        }

        public IEnumerable<Complex> EnumerateScale(double scale)
        {
            var scaleIndex = Array.BinarySearch(_scales, scale);
            if (scaleIndex < 0)
                scaleIndex = ~scaleIndex;
            return EnumerateScale(scaleIndex);
        }

        public IEnumerable<Complex> EnumerateScale(int scaleIndex)
        {
            for (int i = 0; i < _values.GetLength(0); i++)
                yield return _values[i, scaleIndex];
        }

        public IEnumerator<Complex> GetEnumerator()
        {
            for (int i = 0; i < _values.GetLength(0); i++)            
                for (int j = 0; j < _values.GetLength(1); j++)                
                    yield return _values[i, j];    
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
