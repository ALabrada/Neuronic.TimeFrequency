using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics;
using Neuronic.TimeFrequency.Wavelets;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Represents the Continuous Wavelet Transform (CWT).  
    /// </summary>
    public class ContinuousWaveletTransform: IEnumerable<Complex>, IBilinearTimeFrequencyRepresentation
    {
        private readonly Complex[,] _values;
        private readonly double[] _scales;

        private ContinuousWaveletTransform(Complex[,] values, double startTime, double samplingPeriod,
            IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            _values = values;
            SamplingPeriod = samplingPeriod;
            Wavelet = wavelet;
            StartTime = startTime;
            _scales = scales.ToArray();
        }

        /// <summary>
        /// Estimates the CWT for the specified signal using the FFT method.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="scales">The scales.</param>
        /// <returns>The computed CWT.</returns>
        /// <remarks>
        /// This algorithm is based on the one proposed by Dr. Hans-Georg Stark on the book
        /// "Wavelets and Signal Processing: An Application-Based Introduction", Springer (2005).
        /// </remarks>
        public static ContinuousWaveletTransform EstimateUsingFFT(IReadOnlySignal<double> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {            
            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Morlet;

            var yHat = new Complex[signal.Count];
            for (int offset = 0; offset < signal.Count; offset++)
                yHat[offset] = signal[offset];
            yHat.FFT();

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

                psiScale.FFT();

                var factor = 1d / Complex.Sqrt(Math.Abs(scaleArray[scale]));                
                for (int offset = 0; offset < signal.Count; offset++)
                {
                    var trans = Complex.Exp(new Complex(0, -1) * t0 * offset * oms / signal.Count);
                    psiScale[offset] = factor * trans * psiScale[offset] * yHat[offset];
                }

                psiScale.IFFT();

                for (int offset = 0; offset < signal.Count - 1; offset++)
                    values[offset, scale] = psiScale[offset + 1];
            }

            return new ContinuousWaveletTransform(values, signal.Start, signal.SamplingPeriod, wavelet, scaleArray);
        }

        /// <summary>
        /// Estimates the CWT for the specified signal using the FFT method.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="scales">The scales.</param>
        /// <returns>The computed CWT.</returns>
        /// <remarks>
        /// This algorithm is based on the one proposed by Dr. Hans-Georg Stark on the book
        /// "Wavelets and Signal Processing: An Application-Based Introduction", Springer (2005).
        /// </remarks>
        public static ContinuousWaveletTransform EstimateUsingFFT(IReadOnlySignal<float> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            return EstimateUsingFFT(signal.Map(x => (double) x), wavelet, scales);
        }

        /// <summary>
        /// Estimates the CWT for the specified signal using the convolution method.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="scales">The scales.</param>
        /// <returns>The computed CWT.</returns>
        /// <remarks>
        /// This algorithm is based on the <c>cwt</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<double> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            if (wavelet is IWavelet<double> realWavelet)
                return EstimateUsingConvolutions(signal, realWavelet, scales);

            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Morlet;

            var psi = wavelet.EvaluateDomain();
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
                var offset = 1 - (kernel.Count & 1);
                for (int i = 0; i < signal.Count; i++)
                    values[i, scale] = factor * buffer[i + offset];
            }

            return new ContinuousWaveletTransform(values, signal.Start, signal.SamplingPeriod, wavelet, scaleArray);
        }

        /// <summary>
        /// Estimates the CWT for the specified signal using the convolution method.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="scales">The scales.</param>
        /// <returns>The computed CWT.</returns>
        /// <remarks>
        /// This algorithm is based on the <c>cwt</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<float> signal, IWavelet<Complex> wavelet, IEnumerable<double> scales)
        {
            return EstimateUsingConvolutions(signal.Map(x => (double) x), wavelet, scales);
        }

        /// <summary>
        /// Estimates the CWT for the specified signal using the convolution method.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The real-valued wavelet.</param>
        /// <param name="scales">The scales.</param>
        /// <returns>The computed CWT.</returns>
        /// <remarks>
        /// This algorithm is based on the <c>cwt</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<double> signal, IWavelet<double> wavelet, IEnumerable<double> scales)
        {
            scales = scales ?? Enumerable.Range(1, signal.Count).Select(i => signal.SamplingPeriod * i);
            var scaleArray = scales.ToArray();
            wavelet = wavelet ?? Wavelets.Wavelets.Morlet;

            var psi = wavelet.EvaluateDomain();
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
                var offset = 1 - (kernel.Count & 1);
                for (int i = 0; i < signal.Count; i++)
                    values[i, scale] = factor * buffer[i + offset];
            }

            return new ContinuousWaveletTransform(values, signal.Start, signal.SamplingPeriod, (IWavelet<Complex>) wavelet, scaleArray);
        }

        /// <summary>
        /// Estimates the CWT for the specified signal using the convolution method.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The real-valued wavelet.</param>
        /// <param name="scales">The scales.</param>
        /// <returns>The computed CWT.</returns>
        /// <remarks>
        /// This algorithm is based on the <c>cwt</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        public static ContinuousWaveletTransform EstimateUsingConvolutions(IReadOnlySignal<float> signal, IWavelet<double> wavelet, IEnumerable<double> scales)
        {
            return EstimateUsingConvolutions(signal.Map(x => (double) x), wavelet, scales);
        }

        /// <summary>
        /// Gets the offset of the first sample in the time domain.
        /// </summary>
        public double StartTime { get; }

        /// <summary>
        /// Gets the sampling period of the source signal.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the amount of samples in the time domain.
        /// </summary>
        public int SampleCount => _values.GetLength(0);

        /// <summary>
        /// Gets the wavelet function.
        /// </summary>
        public IWavelet<Complex> Wavelet { get; }

        /// <summary>
        /// Gets the estimated scales.
        /// </summary>
        public IList<double> Scales => Array.AsReadOnly(_scales);

        /// <summary>
        /// Gets the amount of samples in the frequency domain.
        /// </summary>
        int IBilinearTimeFrequencyRepresentation.FrequencyCount => _values.GetLength(1);

        /// <summary>
        /// Gets the computed frequencies.
        /// </summary>
        IEnumerable<double> IBilinearTimeFrequencyRepresentation.Frequencies => _scales.Select(s => Wavelet.GetFrequencyOf(s, SamplingPeriod));

        double ITimeFrequencyRepresentation.this[int offset, double frequency]
        {
            get
            {
                var scale = Wavelet.GetScaleFor(frequency, SamplingPeriod);
                var scaleIndex = Array.BinarySearch(_scales, scale);
                if (scaleIndex < 0)
                    scaleIndex = ~scaleIndex;
                return this[offset, scaleIndex].MagnitudeSquared();
            }
        }

        double IBilinearTimeFrequencyRepresentation.this[int offset, int frequencyIndex] => this[offset, frequencyIndex].MagnitudeSquared();

        /// <summary>
        /// Gets the time-frequency content at the specified offset and scale.
        /// </summary>
        /// <param name="offset">The offset.</param>
        /// <param name="scaleIndex">The index of the scale in <see cref="Scales"/>.</param>
        /// <returns>The time-frequency content.</returns>
        public Complex this[int offset, int scaleIndex] => _values[offset, scaleIndex];

        /// <summary>
        /// Gets the time-frequency content at the specified offset and scale.
        /// </summary>
        /// <param name="delay">The time.</param>
        /// <param name="scale">The scale.</param>
        /// <returns>The time-frequency content.</returns>
        public Complex this[double delay, double scale]
        {
            get
            {
                var offset = this.FindClosestOffsetOfTime(delay);
                var scaleIndex = Array.BinarySearch(_scales, scale);
                if (scaleIndex < 0)
                    scaleIndex = ~scaleIndex;
                return this[offset, scaleIndex];
            }
        }

        /// <summary>
        /// Enumerates the computed time domain content for the specified scale.
        /// </summary>
        /// <param name="scale">The scale.</param>
        /// <returns>The time domain content.</returns>
        public IEnumerable<Complex> EnumerateValuesOfScale(double scale)
        {
            var scaleIndex = Array.BinarySearch(_scales, scale);
            if (scaleIndex < 0)
                scaleIndex = ~scaleIndex;
            return EnumerateValuesOfScaleAt(scaleIndex);
        }

        /// <summary>
        /// Enumerates the computed time domain content for the specified scale.
        /// </summary>
        /// <param name="scaleIndex">Index of the scale in <see cref="Scales"/>.</param>
        /// <returns>The time domain content.</returns>
        public IEnumerable<Complex> EnumerateValuesOfScaleAt(int scaleIndex)
        {
            for (int i = 0; i < _values.GetLength(0); i++)
                yield return _values[i, scaleIndex];
        }

        /// <summary>
        /// Returns an enumerator that iterates through the computed time-frequency values.
        /// </summary>
        /// <returns>
        /// An enumerator that can be used to iterate through the collection.
        /// </returns>
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
