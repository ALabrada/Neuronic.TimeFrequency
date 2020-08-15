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
        private readonly Complex[,] _values;
        private readonly double[] _scales;

        private ContinuousWaveletTransform(Complex[,] values, double samplingPeriod, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            _values = values;
            SamplingPeriod = samplingPeriod;
            Wavelet = wavelet;
            _scales = scales.ToArray();
        }

        public static ContinuousWaveletTransform EstimateUsingFFT(IReadOnlySignal<double> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {            
            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Morlet;

            var yHat = new Complex[signal.Count];
            for (int offset = 0; offset < signal.Count; offset++)
                yHat[offset] = signal[offset];
            FourierTransform2.FFT(yHat, FourierTransform.Direction.Forward);

            var t0 = -(signal.Count - 1 - 0.5 * signal.Count) * signal.SamplingPeriod;
            var oms = 2 * Math.PI / signal.SamplingPeriod;
            var values = new Complex[signal.Count - 1, scaleArray.Length];
            var psiScale = new Complex[signal.Count];

            for (int scale = 0; scale < scaleArray.Length; scale++)
            {                
                var psi = new Signal<Complex>(psiScale, t0 / scaleArray[scale], signal.SamplingRate * scaleArray[scale]);
                wavelet.Evaluate(psi);
                psi.Conjugate();
                Array.Reverse(psiScale, 0, psiScale.Length);

                FourierTransform2.FFT(psiScale, FourierTransform.Direction.Forward);

                var factor = 1d / Complex.Sqrt(Math.Abs(scaleArray[scale]));                
                for (int offset = 0; offset < signal.Count; offset++)
                {
                    var trans = Complex.Exp(new Complex(0, -1) * t0 * offset * oms / signal.Count);
                    psiScale[offset] = factor * trans * psiScale[offset] * yHat[offset];
                }

                FourierTransform2.FFT(psiScale, FourierTransform.Direction.Backward);

                for (int offset = 0; offset < signal.Count - 1; offset++)
                    values[offset, scale] = psiScale[offset + 1];
            }

            return new ContinuousWaveletTransform(values, signal.SamplingPeriod, wavelet, scaleArray);
        }

        public static ContinuousWaveletTransform EstimateUsingFFT(IReadOnlySignal<float> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            return EstimateUsingFFT(signal.Map(x => (double) x), wavelet, scales);
        }

        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<double> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            if (wavelet is IWavelet<double> realWavelet)
                return EstimateUsingConvolutions(signal, realWavelet, scales);

            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Morlet;

            var psi = wavelet.Evaluate();
            psi.Integrate();
            psi.Conjugate();

            var values = new Complex[signal.Count, scaleArray.Length];
            var buffer = new Complex[signal.Count + 2];
            var kernel = new List<Complex>(psi.Count);
            for (int scale = 0; scale < scaleArray.Length; scale++)
            {
                var freq = scaleArray[scale] / signal.SamplingPeriod;
                kernel.Clear();
                kernel.AddRange(psi.Sample(freq));
                while (kernel.Count <= 2)
                    kernel.Add(psi[0]);
                kernel.Reverse();

                signal.Convolve(kernel, buffer);
                new Signal<Complex>(buffer).Differentiate();

                var factor = -Math.Sqrt(scaleArray[scale]);
                for (int i = 0; i < signal.Count; i++)                
                    values[i, scale] = factor * buffer[i + 1];                
            }

            return new ContinuousWaveletTransform(values, signal.SamplingPeriod, wavelet, scaleArray);
        }

        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<float> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            return EstimateUsingConvolutions(signal.Map(x => (double) x), wavelet, scales);
        }

        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<double> signal, IWavelet<double> wavelet, IEnumerable<double> scales)
        {
            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Morlet;

            var psi = wavelet.Evaluate();
            psi.Integrate();

            var values = new Complex[signal.Count, scaleArray.Length];
            var buffer = new double[signal.Count + 2];
            var kernel = new List<double>(psi.Count);
            for (int scale = 0; scale < scaleArray.Length; scale++)
            {
                var freq = scaleArray[scale] / signal.SamplingPeriod;
                kernel.Clear();
                kernel.AddRange(psi.Sample(freq));
                while (kernel.Count <= 2)
                    kernel.Add(psi[0]);
                kernel.Reverse();

                signal.Convolve(kernel, buffer);
                new Signal<double>(buffer).Differentiate();

                var factor = -Math.Sqrt(scaleArray[scale]);
                for (int i = 0; i < signal.Count; i++)
                    values[i, scale] = factor * buffer[i + 1];
            }

            return new ContinuousWaveletTransform(values, signal.SamplingPeriod, (IWavelet<Complex>) wavelet, scaleArray);
        }

        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<float> signal, IWavelet<double> wavelet, IEnumerable<double> scales)
        {
            return EstimateUsingConvolutions(signal.Map(x => (double) x), wavelet, scales);
        }

        public double SamplingPeriod { get; }

        public IWavelet<Complex> Wavelet { get; }

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
