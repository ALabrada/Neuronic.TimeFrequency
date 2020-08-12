using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;
using Neuronic.TimeFrequency.Wavelets;

namespace Neuronic.TimeFrequency
{
    public class ContinuousWaveletTransform: IEnumerable<Complex>
    {
        private Complex[,] _values;
        private double[] _scales;

        private ContinuousWaveletTransform(Complex[,] values, double samplingPeriod, IWavelet wavelet, IEnumerable<double> scales)
        {
            _values = values;
            SamplingPeriod = samplingPeriod;
            Wavelet = wavelet;
            _scales = scales.ToArray();
        }

        public static ContinuousWaveletTransform EstimateUsingFFT(Signal<float> signal, IWavelet wavelet, IEnumerable<double> scales)
        {            
            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Haar;

            var general = new Complex[signal.Count];
            for (int offset = 0; offset < signal.Count; offset++)
                general[offset] = signal.Count;
            FourierTransform2.DFT(general, FourierTransform.Direction.Forward);


            var values = new Complex[signal.Count, scaleArray.Length];
            var scaled = new Complex[signal.Count];

            for (int scale = 0; scale < scaleArray.Length; scale++)
            {                
                var dt = (signal.Count - 1) * signal.SamplingPeriod / scaleArray[scale];
                var phi = new Signal<Complex>(scaled, -dt, signal.SamplingRate);
                wavelet.Evaluate(phi);
                Array.Reverse(scaled, 0, scaled.Length);

                for (int offset = 0; offset < signal.Count; offset++)
                    scaled[offset] = Complex.Conjugate(scaled[offset]);

                FourierTransform2.DFT(scaled, FourierTransform.Direction.Forward);

                var factor = 1d / Complex.Sqrt(wavelet.Energy * Math.Abs(scaleArray[scale]));
                for (int offset = 0; offset < signal.Count; offset++)
                    scaled[offset] = factor * scaled[offset] * general[offset];

                FourierTransform2.DFT(scaled, FourierTransform.Direction.Backward);

                for (int offset = 0; offset < signal.Count; offset++)
                    values[offset, scale] = scaled[offset];
            }

            return new ContinuousWaveletTransform(values, signal.SamplingPeriod, wavelet, scaleArray);
        }

        public static ContinuousWaveletTransform EstimateAsMatlab(Signal<float> signal, IWavelet wavelet, IEnumerable<double> scales)
        {
            throw new NotImplementedException();
        }

        public double SamplingPeriod { get; }

        public IWavelet Wavelet { get; }

#if !NET40
        public IReadOnlyList<double> Scales => _scales;
#else
        public IList<double> Scales => _scales;
#endif

        public Complex this[int offset, int scaleIndex] => _values[offset, scaleIndex];

        public Complex this[double delay, double scale]
        {
            get
            {
                var offset = (int)Math.Round(delay / SamplingPeriod);
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
