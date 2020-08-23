using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Base class for Lag Independent kernels.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Kernels.DopplerLagKernel" />
    public abstract class DopplerKernel : DopplerLagKernel
    {
        /// <summary>
        /// Gets or sets the window function samples.
        /// </summary>
        public double[] Window { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether to use time or Doppler domain.
        /// </summary>
        public bool UseDopplerDomain { get; set; }

        /// <summary>
        /// Evaluates lag kernel.
        /// </summary>
        /// <param name="n">The number of samples in the Doppler domain.</param>
        /// <param name="nTime">The number of samples in the time domain.</param>
        /// <param name="output">The output buffer.</param>
        /// <param name="start">The start index in <paramref name="output"/>.</param>
        /// <returns>The number of samples.</returns>
        protected virtual int Evaluate(int n, int nTime, double[] output, int start)
        {
            if (UseDopplerDomain)
                return EvaluateInDopplerDomain(nTime, output, start);
            return EvaluateInTimeDomain(n, output, start);
        }

        /// <summary>
        /// Evaluates lag kernel.
        /// </summary>
        /// <param name="n">The number of samples in the Doppler domain.</param>
        /// <param name="nTime">The number of samples in the time domain.</param>
        /// <returns>The lag kernel.</returns>
        protected virtual IList<double> Evaluate(int n, int nTime)
        {
            var output = new double[Math.Max(n, nTime)];
            var count = Evaluate(n, nTime, output, 0);
            return new ArraySegment<double>(output, 0, count);
        }

        private int EvaluateInDopplerDomain(int nTime, double[] output, int start)
        {
            var win = Window ?? throw new InvalidOperationException("No window selected.");
            var q = win.Length;

            if (q > nTime / 2)
                throw new ArgumentOutOfRangeException(nameof(nTime));

            var padWin = new double[nTime / 2];
            ShiftAndPad(win, padWin, 0, padWin.Length);
            if (nTime / 2 > q && q % 2 == 0)
                q++;

            var g = output;
            for (int i = 0; i <= q / 2; i++)
                g[start + i] = padWin[i];
            for (int i = 1; i <= q / 2 + (q & 1) - 1; i++)
                g[start + q - i] = padWin[nTime / 2 - i];
            return q;
        }

        private int EvaluateInTimeDomain(int n, double[] output, int start)
        {
            var win = Window ?? throw new InvalidOperationException("No window selected.");
            var re = output;
            ShiftAndPad(win, output, start, n);

            var im = new double[n];
            FourierTransform2.FFT(re, im, FourierTransform.Direction.Forward);

            var win0 = new Complex(re[0], im[0]);
            re[0] = 1;
            for (int i = 1; i < re.Length; i++)
            {
                var c = new Complex(re[i], im[i]);
                re[i] = (c / win0).Real;
            }

            return n;
        }

        private static void ShiftAndPad(double[] win, double[] padWin, int start, int count)
        {
            var q = win.Length;
            // Shift (q+1)/2
            Array.Copy(win, 0, padWin, start + count - q / 2, q / 2);
            Array.Copy(win, q / 2, padWin, start, (q + 1) / 2);

            if (win.Length % 2 == 0)
                padWin[start + count - q / 2] = padWin[start + q / 2] = win[q / 2] / 2d;
        }

        private double[] PadWindow(double[] win, int length)
        {
            if (length < win.Length) throw new ArgumentOutOfRangeException(nameof(length));
            if (length == win.Length) return win;
            var result = new double[length];
            PadWindow(win, result);
            return result;
        }

        /// <summary>
        /// Pads the window to the specified length.
        /// </summary>
        /// <param name="src">The source window.</param>
        /// <param name="dst">The destination window.</param>
        /// <exception cref="ArgumentException">The source does not fit in the destination. - src</exception>
        protected void PadWindow(IList<double> src, IList<double> dst)
        {
            if (dst.Count < src.Count) throw new ArgumentException("The source does not fit in the destination.", nameof(src));

            if (src.Count % 2 == 1)
            {
                for (int i = 0; i <= src.Count / 2; i++)
                    dst[i] = src[i];
                for (int i = 1; i <= src.Count / 2; i++)
                    dst[dst.Count - i] = src[src.Count - i];
            }
            else
            {
                for (int i = 0; i < src.Count / 2; i++)
                    dst[i] = src[i];
                dst[src.Count / 2] = src[src.Count / 2] / 2d;

                for (int i = 1; i < src.Count / 2; i++)
                    dst[dst.Count - i] = src[src.Count - i];
                dst[dst.Count - src.Count / 2] = src[src.Count / 2] / 2d;
            }


            for (int i = src.Count / 2 + 1; i < dst.Count - src.Count / 2; i++)
                dst[i] = 0;
        }
    }
}