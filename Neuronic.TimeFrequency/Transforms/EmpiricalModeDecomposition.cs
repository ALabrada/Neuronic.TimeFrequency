using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Empirical Mode Decomposition (EMD) in Intrinsic Mode Functions (IMF).
    /// </summary>
    public class EmpiricalModeDecomposition: IReadOnlyList<IReadOnlySignal<double>>
    {
        private readonly IList<double[]> _imfs;

        private EmpiricalModeDecomposition(IList<double[]> imfs, double samplingPeriod)
        {
            _imfs = imfs;
            SamplingPeriod = samplingPeriod;
        }

        private static void GetLocalExtremeValues(IReadOnlyList<double> signal, List<DoublePoint> max, List<DoublePoint> min)
        {
            max.Clear();
            min.Clear();
            if (signal.Count <= 2)
                return;

            var bounds = signal.Skip(1).Aggregate(new {min = signal[0], max = signal[0]}, (b, x) => new 
            {
                min = Math.Min(b.min, x), max = Math.Max(b.max, x)
            });
            var maxGap = bounds.max - bounds.min;

            void Interpolate(int i, IList<DoublePoint> output)
            {
                var ya = signal[i - 1];
                var yb = signal[i];
                var yc = signal[i + 1];
                DoublePoint point;

                var discriminat = -4 * yb + 2 * ya + 2 * yc;
                point.X = discriminat == 0 ? i : i + (ya - yc) / discriminat;

                discriminat = -16 * yb + 8 * ya + 8 * yc;
                point.Y = discriminat == 0 ? yb : yb + (2 * yc * ya - ya * ya - yc * yc) / discriminat;

                var previous = output.Count == 0
                    ? new DoublePoint(Double.NegativeInfinity, Double.NegativeInfinity)
                    : output[output.Count - 1];
                if (point.X == previous.X || Math.Abs(point.Y - previous.Y) / Math.Abs(point.X - previous.Y) > 2 * maxGap)
                {
                    point.X = (point.X + previous.X) / 2;
                    point.Y = Math.Max(point.Y, previous.Y);
                    output[output.Count - 1] = point;
                }
                else
                    output.Add(point);
            }

            for (int i = 1; i < signal.Count - 1; i++)
            {
                var ya = signal[i - 1];
                var yb = signal[i];
                var yc = signal[i + 1];

                if (yb > yc && yb >= ya || yb >= yc && yb > ya)
                    Interpolate(i, max);
                else if (yb < yc && yb <= ya || yb <= yc && yb < ya)
                    Interpolate(i, min);
            }

            if (max.Count > 0)
            {
                if (signal[0] >= max[0].Y)
                    max.Insert(0, new DoublePoint(0, signal[0]));
                if (signal[signal.Count - 1] >= max[max.Count - 1].Y) 
                    max.Add(new DoublePoint(signal.Count - 1, signal[signal.Count - 1]));
            }
            else
            {
                if (signal[0] > signal[1])
                    max.Insert(0, new DoublePoint(0, signal[0]));
                if (signal[signal.Count - 1] > signal[signal.Count - 2])
                    max.Add(new DoublePoint(signal.Count - 1, signal[signal.Count - 1]));
            }

            if (min.Count > 0)
            {
                if (signal[0] <= min[0].Y)
                    min.Insert(0, new DoublePoint(0, signal[0]));
                if (signal[signal.Count - 1] <= min[min.Count - 1].Y)
                    min.Add(new DoublePoint(signal.Count - 1, signal[signal.Count - 1]));
            }
            else
            {
                if (signal[0] < signal[1])
                    min.Insert(0, new DoublePoint(0, signal[0]));
                if (signal[signal.Count - 1] < signal[signal.Count - 2])
                    min.Add(new DoublePoint(signal.Count - 1, signal[signal.Count - 1]));
            }

            max.Insert(0, new DoublePoint(-min[0].X, max[0].Y));
            max.Add(new DoublePoint(2 * (signal.Count - 1) - min[min.Count - 1].X, max[max.Count - 1].Y));

            min.Insert(0, new DoublePoint(-max[1].X, min[0].Y));
            min.Add(new DoublePoint(2 * (signal.Count - 1) - max[max.Count - 2].X, min[min.Count - 1].Y));

            max.Sort((p1, p2) => p1.X.CompareTo(p2.X));
            min.Sort((p1, p2) => p1.X.CompareTo(p2.X));
        }

        private static IEnumerable<double> SplineInterpolation(IReadOnlyList<DoublePoint> source, IEnumerable<double> eval)
        {
            Func<double, double> interpolation;
            var n = source.Count;
            if (n < 2)
                throw new ArgumentException("2 points at least are required.", nameof(source));

            var xs = new double[n];
            var ys = new double[n];

            source.Select(v => v.X + 1).CopyTo(xs, 0);
            source.Select(v => v.Y).CopyTo(ys, 0);

            double dx(int i) => xs[i] - xs[i - 1];
            double divdif(int i) => (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);

            if (n == 2)
            {
                //var temp = divdif(1);
                //ys[1] = ys[0];
                //ys[0] = temp;

                interpolation = StepInterpolation.InterpolateSorted(xs, ys).Interpolate;
            }
            else if (n == 3)
            {
                ys[2] = divdif(2);
                ys[1] = divdif(1);
                ys[2] = (ys[2] - ys[1]) / (xs[2] - xs[0]);
                ys[1] -= ys[2] * dx(1);

                var p = new Polynomial(ys);
                interpolation = p.Evaluate;
            }
            else
            {
                var b = new DenseVector(n);
                for (int i = 1; i < n - 1; i++)
                    b[i] = 3d * (dx(i + 1) * divdif(i) + dx(i) * divdif(i + 1));

                var x31 = xs[2] - xs[0];
                var xn = xs[xs.Length - 1] - xs[xs.Length - 3];
                b[0] = ((dx(1) + 2 * x31) * dx(2) * divdif(1) + dx(1) * dx(1) * divdif(2)) / x31;
                b[n - 1] = (dx(n - 1) * dx(n - 1) * divdif(n - 2) + (2 * xn + dx(n - 1)) * dx(n - 2) * divdif(n - 1)) / xn;

                var c = new SparseMatrix(n, n)
                {
                    [0, 1] = x31,
                    [0, 0] = dx(2),
                    [n - 1, n - 1] = dx(n - 2),
                    [n - 1, n - 2] = xn
                };
                for (int i = 1; i <= n - 2; i++)
                {
                    c[i, i + 1] = dx(i);
                    c[i, i] = 2d * (dx(i + 1) + dx(i));
                    c[i, i - 1] = dx(i + 1);
                }

                var slopes = c.Solve(b).ToArray();

                var spline = CubicSpline.InterpolateHermiteSorted(xs, ys, slopes);
                interpolation = spline.Interpolate;
            }
            
            foreach (var x in eval)
                yield return interpolation(x + 1);
        }

        /// <summary>
        /// Estimates the EMD of the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="outerStop">The stop criteria of the outer loop. Default is <see cref="ResidualStopCriteria"/>.</param>
        /// <param name="innerStop">The stop criteria of the inner loop, responsible for obtaining each IMF. Default is <see cref="ResolutionStopCriteria"/>.</param>
        /// <param name="alpha">Gradient step size, in the range (0, 1].</param>
        /// <returns>The EMD of <paramref name="signal"/>.</returns>
        /// <remarks>
        /// This algorithm is based on the work of Rato, R. T., Ortigueira, M. D. and Batista, A. G., published in the article
        /// "On the HHT, its problems, and some solutions",
        /// Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.
        /// </remarks>
        public static EmpiricalModeDecomposition Estimate(IReadOnlySignal<double> signal, 
            IStopCriteria<OuterState> outerStop = null, IStopCriteria<InnerState> innerStop = null, 
            double alpha = 1d)
        {
            if (alpha <= 0 || alpha > 1) throw new ArgumentOutOfRangeException(nameof(alpha));

            var localMax = new List<DoublePoint>(signal.Count);
            var localMin = new List<DoublePoint>(signal.Count);
            var results = new List<double[]>();

            var samples = new Signal<double>(new double[signal.Count], signal.Start, signal.SamplingRate);
            signal.CopyTo(samples.Samples, 0);
            
            var bias = new Signal<double>(new double[signal.Count], signal.Start, signal.SamplingRate);

            var energy = signal.Sum(x => x * x);

            outerStop = outerStop ?? new ResidualStopCriteria(energy);
            innerStop = innerStop ?? new ResolutionStopCriteria();

            for (int itOut = 0; !outerStop.ShouldStop(new OuterState(itOut, samples, results)); itOut++)
            {
                var imf = new Signal<double>(new double[signal.Count], signal.Start, signal.SamplingRate);
                samples.CopyTo(imf.Samples, 0);

                GetLocalExtremeValues(imf, localMax, localMin);

                Debug.Assert(Math.Abs(localMax.Count - localMin.Count) < 2, "Max-Min count mismatch.");

                var top = SplineInterpolation(localMax, Enumerable.Range(0, signal.Count).Select(x => (double) x));
                var bottom = SplineInterpolation(localMin, Enumerable.Range(0, signal.Count).Select(x => (double) x));
                top.Zip(bottom, (t, b) => (t + b) / 2).CopyTo(bias.Samples, 0);

                for (int itIn = 0; !innerStop.ShouldStop(new InnerState(itIn, imf, bias)); itIn++)
                {
                    for (int i = 0; i < imf.Count; i++)
                        imf[i] -= alpha * bias[i];

                    GetLocalExtremeValues(imf, localMax, localMin);

                    top = SplineInterpolation(localMax, Enumerable.Range(0, signal.Count).Select(x => (double)x));
                    bottom = SplineInterpolation(localMin, Enumerable.Range(0, signal.Count).Select(x => (double)x));
                    top.Zip(bottom, (t, b) => (t + b) / 2).CopyTo(bias.Samples, 0);
                }

                results.Add(imf.Samples);

                for (int i = 0; i < imf.Count; i++)
                    samples[i] -= imf[i];
            }

            if (samples.Sum(x => x * x) / energy > 10e-12)
                results.Add(samples.Samples);

            return new EmpiricalModeDecomposition(results, signal.SamplingPeriod);
        }

        /// <summary>
        /// Estimates the EMD of the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="outerStop">The stop criteria of the outer loop.</param>
        /// <param name="innerStop">The stop criteria of the inner loop, responsible for obtaining each IMF.</param>
        /// <param name="alpha">Gradient step size, in the range (0, 1].</param>
        /// <returns>The EMD of <paramref name="signal"/>.</returns>
        /// <remarks>
        /// This algorithm is based on the work of Rato, R. T., Ortigueira, M. D. and Batista, A. G., published in the article
        /// "On the HHT, its problems, and some solutions",
        /// Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.
        /// </remarks>
        public static EmpiricalModeDecomposition Estimate(IReadOnlySignal<float> signal,
            IStopCriteria<OuterState> outerStop = null, IStopCriteria<InnerState> innerStop = null,
            double alpha = 1d)
        {
            return Estimate(signal.Map(x => (double) x), outerStop, innerStop, alpha);
        }

        /// <summary>
        /// Gets the sampling period.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the amount of IMFs.
        /// </summary>
        public int Count => _imfs.Count;

        /// <summary>
        /// Gets the IMF at the specified index.
        /// </summary>
        /// <param name="index">The index.</param>
        /// <returns>The IMF.</returns>
        public IReadOnlySignal<double> this[int index] => new Signal<double>(_imfs[index], fs: 1/SamplingPeriod);

        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>
        /// An enumerator that can be used to iterate through the collection.
        /// </returns>
        public IEnumerator<IReadOnlySignal<double>> GetEnumerator()
        {
            return Enumerable.Range(0, Count).Select(i => this[i]).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        /// <summary>
        /// Performs the spectral analysis of the EMD using local AR approximations.
        /// </summary>
        /// <param name="windowLength">Length of the analysis window.</param>
        /// <returns>The spectral analysis of the EMD.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="windowLength"/> is not positive.</exception>
        /// <remarks>
        /// This algorithm is based on the method proposed by Rato, R. T., Ortigueira, M. D. and Batista, A. G. in the article
        /// "On the HHT, its problems, and some solutions",
        /// Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.
        /// </remarks>
        public SpectralAnalysis LocalSpectralAnalysis(int windowLength = 40)
        {
            var n = _imfs[0].Length;
            if (windowLength <= 0 || windowLength > n) throw new ArgumentOutOfRangeException(nameof(windowLength));
            var max = new List<DoublePoint>(n);
            var min = new List<DoublePoint>(n);
            var components = new List<SpectralAnalysis.MonocomponentSignal>(_imfs.Count);

            foreach (var imf in _imfs)
            {
                GetLocalExtremeValues(imf, max, min);
                
                var amplitude = new double[n];
                SplineInterpolation(max, Enumerable.Range(0, imf.Length).Select(x => (double)x)).CopyTo(amplitude, 0);

                double FM(int index) => imf[index] / amplitude[index];

                var frequency = new double[n - windowLength + 1];

                for (int i = 0; i < frequency.Length; i++)
                {
                    var sqSum = 0d;
                    var sum = 0d;
                    for (int j = 1, k = i + j; j < windowLength - 1; j++, k++)
                    {
                        var x = FM(k);
                        sum += x * (FM(k - 1) + FM(k + 1));
                        sqSum += 2 * x * x;
                    }

                    frequency[i] = Math.Acos(sum / sqSum) / (2 * Math.PI * SamplingPeriod);
                }

                components.Add(new SpectralAnalysis.MonocomponentSignal(amplitude, frequency));
            }

            return new SpectralAnalysis(components, SamplingPeriod, n);
        }

        /// <summary>
        /// Performs the Hilbert Spectral Analysis of the EMD.
        /// </summary>
        /// <returns>The Hilbert-Huang Transform of the signal.</returns>
        /// <remarks>
        /// This algorithm is based on Matlab.
        /// </remarks>
        public SpectralAnalysis HilbertSpectralAnalysis()
        {
            var fs = 1d / SamplingPeriod;
            var ws = fs / (2 * Math.PI);

            var components = new List<SpectralAnalysis.MonocomponentSignal>(_imfs.Count);
            var z = new Signal<Complex>(new Complex[_imfs[0].Length]);

            foreach (var imf in _imfs)
            {
                imf.Select(x => (Complex) x).CopyTo(z.Samples, 0);

                z.HilbertTransform();

                var amplitude = new double[z.Count];
                z.Select(c => c.Magnitude).CopyTo(amplitude, 0);

                var frequency = new Signal<double>(new double[z.Count], fs: fs);
                z.Select(c => c.Phase).CopyTo(frequency.Samples, 0);

                frequency.Unwrap();
                frequency.Differentiate();

                for (int i = 0; i < frequency.Count; i++)
                    frequency[i] *= ws;

                var component = new SpectralAnalysis.MonocomponentSignal(amplitude, 
                    new ArraySegment<double>(frequency.Samples, 0, frequency.Count - 1));
                components.Add(component);
            }

            return new SpectralAnalysis(components, SamplingPeriod, z.Count);
        }

        private struct DoublePoint
        {
            public DoublePoint(double x, double y)
            {
                X = x;
                Y = y;
            }

            public double X;
            public double Y;
        }

        /// <summary>
        /// Contains the state of the outer loop of the iterative algorithm used in
        /// <see cref="Estimate(IReadOnlySignal{double},IStopCriteria{OuterState},IStopCriteria{InnerState},double)"/>.
        /// </summary>
        public struct OuterState
        {
            internal OuterState(int iteration, IReadOnlySignal<double> signal, IReadOnlyList<IReadOnlyList<double>> imfs)
            {
                Signal = signal;
                IntrinsicModeFunctions = imfs;
                Iteration = iteration;
            }

            /// <summary>
            /// Gets the signal.
            /// </summary>
            public IReadOnlySignal<double> Signal { get; }

            /// <summary>
            /// Gets the obtained intrinsic mode functions.
            /// </summary>
            public IReadOnlyList<IReadOnlyList<double>> IntrinsicModeFunctions { get; }

            /// <summary>
            /// Gets the iteration.
            /// </summary>
            public int Iteration { get; }
        }

        /// <summary>
        /// Contains the state of the inner loop of the iterative algorithm used in
        /// <see cref="Estimate(IReadOnlySignal{double},IStopCriteria{OuterState},IStopCriteria{InnerState},double)"/>.
        /// </summary>
        public struct InnerState
        {
            internal InnerState(int iteration, IReadOnlySignal<double> imf, IReadOnlySignal<double> bias)
            {
                IntrinsicModeFunction = imf;
                Bias = bias;
                Iteration = iteration;
            }

            /// <summary>
            /// Gets the obtained intrinsic mode function.
            /// </summary>
            public IReadOnlySignal<double> IntrinsicModeFunction { get; }

            /// <summary>
            /// Gets the average of the signal envelopes.
            /// </summary>
            public IReadOnlySignal<double> Bias { get; }

            /// <summary>
            /// Gets the iteration.
            /// </summary>
            public int Iteration { get; }
        }
    }

    /// <summary>
    /// The default stop criteria for the outer loop of
    /// <see cref="EmpiricalModeDecomposition.Estimate(IReadOnlySignal{double},IStopCriteria{EmpiricalModeDecomposition.OuterState},IStopCriteria{EmpiricalModeDecomposition.InnerState},double)"/>
    /// </summary>
    public class ResidualStopCriteria: IStopCriteria<EmpiricalModeDecomposition.OuterState>
    {
        private readonly double _energy;

        /// <summary>
        /// Initializes a new instance of the <see cref="ResidualStopCriteria"/> class.
        /// </summary>
        /// <param name="originalSignal">The original signal.</param>
        /// <param name="residual">The residual energy limit in dB. Normally between 40 and 60 dB.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when <paramref name="residual"/> is not a positive value.
        /// </exception>
        public ResidualStopCriteria(IReadOnlySignal<double> originalSignal, double residual = 50d)
        : this (originalSignal.Sum(x => x*x), residual)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ResidualStopCriteria"/> class.
        /// </summary>
        /// <param name="originalSignal">The original signal.</param>
        /// <param name="residual">The residual energy limit in dB. Normally between 40 and 60 dB.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when <paramref name="residual"/> is not a positive value.
        /// </exception>
        public ResidualStopCriteria(IReadOnlySignal<float> originalSignal, double residual = 50d)
            : this (originalSignal.Map(x => (double) x), residual)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ResidualStopCriteria"/> class.
        /// </summary>
        /// <param name="energy">The energy.</param>
        /// <param name="residual">The residual energy in dB. Normally between 40 and 60 dB.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when <paramref name="residual"/> is not a positive value.
        /// </exception>
        public ResidualStopCriteria(double energy, double residual = 50d)
        {
            if (residual <= 0) throw new ArgumentOutOfRangeException(nameof(residual));
            if (energy < 0) throw new ArgumentOutOfRangeException(nameof(energy));
            Residual = residual;
            _energy = energy;
        }

        /// <summary>
        /// Gets the residual energy limit.
        /// </summary>
        public double Residual { get; }

        /// <summary>
        /// Determines if the algorithm should stop iterating.
        /// </summary>
        /// <param name="state">The current state.</param>
        /// <returns>
        ///   <c>true</c> if the algorithm should stop; otherwise, <c>false</c>.
        /// </returns>
        bool IStopCriteria<EmpiricalModeDecomposition.OuterState>.ShouldStop(EmpiricalModeDecomposition.OuterState state)
        {
            var oscCount = state.Signal.CountOscilations();
            if (oscCount <= 2)
                return true;

            var sampleEnergy = state.Signal.Sum(x => x * x);
            var currentResidual = sampleEnergy > 0 ? 10 * Math.Log10(_energy / sampleEnergy) : double.PositiveInfinity;

            return currentResidual >= Residual;
        }
    }


    /// <summary>
    /// The default stop criteria for the inner loop of
    /// <see cref="EmpiricalModeDecomposition.Estimate(IReadOnlySignal{double},IStopCriteria{EmpiricalModeDecomposition.OuterState},IStopCriteria{EmpiricalModeDecomposition.InnerState},double)"/>
    /// </summary>
    public class ResolutionStopCriteria : IStopCriteria<EmpiricalModeDecomposition.InnerState>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="ResolutionStopCriteria"/> class.
        /// </summary>
        /// <param name="resolution">The resolution limit in dB. Normally between 40 and 60 dB.</param>
        /// <exception cref="ArgumentOutOfRangeException">resolution</exception>
        public ResolutionStopCriteria(double resolution = 50d)
        {
            if (resolution <= 0) throw new ArgumentOutOfRangeException(nameof(resolution));
            Resolution = resolution;
        }

        /// <summary>
        /// Gets the resolution limit.
        /// </summary>
        public double Resolution { get; }

        /// <summary>
        /// Determines if the algorithm should stop iterating.
        /// </summary>
        /// <param name="state">The current state.</param>
        /// <returns>
        ///   <c>true</c> if the algorithm should stop; otherwise, <c>false</c>.
        /// </returns>
        bool IStopCriteria<EmpiricalModeDecomposition.InnerState>.ShouldStop(EmpiricalModeDecomposition.InnerState state)
        {
            var currentEnergy = state.IntrinsicModeFunction.Sum(x => x * x);
            var biasEnergy = state.Bias.Sum(x => x * x);
            var currentResolution = biasEnergy > 0
                ? 10 * Math.Log10(currentEnergy / biasEnergy)
                : double.PositiveInfinity;
            return currentResolution > Resolution;
        }
    }

    /// <summary>
    /// A stop criteria for 
    /// <see cref="EmpiricalModeDecomposition.Estimate(IReadOnlySignal{double},IStopCriteria{EmpiricalModeDecomposition.OuterState},IStopCriteria{EmpiricalModeDecomposition.InnerState},double)"/>
    /// with a fixed number of iterations.
    /// </summary>
    public class FixedStopCriteria : 
        IStopCriteria<EmpiricalModeDecomposition.InnerState>,
        IStopCriteria<EmpiricalModeDecomposition.OuterState>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="FixedStopCriteria"/> class.
        /// </summary>
        /// <param name="iterationCount">The iteration count.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="iterationCount"/> is not positive.</exception>
        public FixedStopCriteria(int iterationCount)
        {
            if (iterationCount <= 0) throw new ArgumentOutOfRangeException(nameof(iterationCount));
            IterationCount = iterationCount;
        }

        /// <summary>
        /// Gets the iteration count.
        /// </summary>
        public int IterationCount { get; }

        bool IStopCriteria<EmpiricalModeDecomposition.InnerState>.ShouldStop(EmpiricalModeDecomposition.InnerState state)
        {
            return state.Iteration >= IterationCount;
        }

        bool IStopCriteria<EmpiricalModeDecomposition.OuterState>.ShouldStop(EmpiricalModeDecomposition.OuterState state)
        {
            return state.Iteration >= IterationCount;
        }
    }
}