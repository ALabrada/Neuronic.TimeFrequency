using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Accord;
using Accord.Diagnostics;
using MathNet.Numerics.Interpolation;

namespace Neuronic.TimeFrequency
{
    public class EmpiricalModeDecomposition: IReadOnlyList<IReadOnlySignal<double>>
    {
        private readonly IList<double[]> _imfs;

        public EmpiricalModeDecomposition(IList<double[]> imfs, double samplingPeriod)
        {
            _imfs = imfs;
            SamplingPeriod = samplingPeriod;
        }

        private static void GetLocalExtremeValues(IReadOnlyList<double> signal, IList<DoublePoint> max, IList<DoublePoint> min)
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
                    output[max.Count - 1] = point;
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
        }

        private static IEnumerable<double> SplineInterpolation(IReadOnlyList<DoublePoint> source, IEnumerable<double> eval)
        {
            var spline = CubicSpline.InterpolateNatural(source.Select(v => v.X), source.Select(v => v.Y));
            foreach (var x in eval)
                yield return spline.Interpolate(x);
        }

        public static EmpiricalModeDecomposition Estimate(IReadOnlySignal<double> signal, 
            IStopCriteria<OuterState> outerStop = null, IStopCriteria<InnerState> innerStop = null, 
            double alpha = 1d)
        {
            if (alpha <= 0 || alpha > 1) throw new ArgumentOutOfRangeException(nameof(alpha));
            outerStop = outerStop ?? new ResidualStopCriteria(signal);
            innerStop = innerStop ?? new ResolutionStopCriteria();

            var localMax = new List<DoublePoint>(signal.Count);
            var localMin = new List<DoublePoint>(signal.Count);
            var results = new List<double[]>();

            var samples = new Signal<double>(new double[signal.Count], signal.Start, signal.SamplingRate);
            signal.CopyTo(samples.Samples, 0);


            var imf = new Signal<double>(new double[signal.Count], signal.Start, signal.SamplingRate);
            var bias = new Signal<double>(new double[signal.Count], signal.Start, signal.SamplingRate);

            var energy = signal.Sum(x => x * x);
            for (int itOut = 0; !outerStop.ShouldStop(new OuterState(itOut, samples)); itOut++)
            {
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

        public static EmpiricalModeDecomposition Estimate(IReadOnlySignal<float> signal,
            IStopCriteria<OuterState> outerStop = null, IStopCriteria<InnerState> innerStop = null,
            double alpha = 1d)
        {
            return Estimate(signal.Map(x => (double) x), outerStop, innerStop, alpha);
        }

        public double SamplingPeriod { get; }

        public int Count => _imfs.Count;

        public IReadOnlySignal<double> this[int index] => new Signal<double>(_imfs[index], fs: 1/SamplingPeriod);

        public IEnumerator<IReadOnlySignal<double>> GetEnumerator()
        {
            return Enumerable.Range(0, Count).Select(i => this[i]).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public struct OuterState
        {
            public OuterState(int iteration, IReadOnlySignal<double> imf)
            {
                IntrinsicModeFunction = imf;
                Iteration = iteration;
            }

            public IReadOnlySignal<double> IntrinsicModeFunction { get; }

            public int Iteration { get; }
        }

        public struct InnerState
        {
            internal InnerState(int iteration, IReadOnlySignal<double> imf, IReadOnlySignal<double> bias)
            {
                IntrinsicModeFunction = imf;
                Bias = bias;
                Iteration = iteration;
            }

            public IReadOnlySignal<double> IntrinsicModeFunction { get; }

            public IReadOnlySignal<double> Bias { get; }

            public int Iteration { get; }
        }
    }

    public class ResidualStopCriteria: IStopCriteria<EmpiricalModeDecomposition.OuterState>
    {
        private readonly double _energy;

        public ResidualStopCriteria(IReadOnlySignal<double> originalSignal, double residual = 50d)
        {
            if (originalSignal == null) throw new ArgumentNullException(nameof(originalSignal));
            if (residual <= 0) throw new ArgumentOutOfRangeException(nameof(residual));
            Residual = residual;
            _energy = originalSignal.Sum(x => x * x);
        }

        public ResidualStopCriteria(IReadOnlySignal<float> originalSignal, double residual = 50d)
            : this (originalSignal.Map(x => (double) x), residual)
        {
        }

        public double Residual { get; }

        public bool ShouldStop(EmpiricalModeDecomposition.OuterState state)
        {
            var oscCount = state.IntrinsicModeFunction.CountOscilations();
            if (oscCount <= 2)
                return true;

            var sampleEnergy = state.IntrinsicModeFunction.Sum(x => x * x);
            var currentResidual = sampleEnergy > 0 ? 10 * Math.Log10(_energy / sampleEnergy) : double.PositiveInfinity;

            return currentResidual >= Residual;
        }
    }

    public class ResolutionStopCriteria : IStopCriteria<EmpiricalModeDecomposition.InnerState>
    {
        public ResolutionStopCriteria(double resolution = 50d)
        {
            if (resolution <= 0) throw new ArgumentOutOfRangeException(nameof(resolution));
            Resolution = resolution;
        }

        public double Resolution { get; }

        public bool ShouldStop(EmpiricalModeDecomposition.InnerState state)
        {
            var currentEnergy = state.IntrinsicModeFunction.Sum(x => x * x);
            var biasEnergy = state.Bias.Sum(x => x * x);
            var currentResolution = biasEnergy > 0
                ? 10 * Math.Log10(currentEnergy / biasEnergy)
                : double.PositiveInfinity;
            return currentResolution > Resolution;
        }
    }

    public class FixedStopCriteria : 
        IStopCriteria<EmpiricalModeDecomposition.InnerState>,
        IStopCriteria<EmpiricalModeDecomposition.OuterState>
    {
        public FixedStopCriteria(int iterationCount)
        {
            if (iterationCount <= 0) throw new ArgumentOutOfRangeException(nameof(iterationCount));
            IterationCount = iterationCount;
        }

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