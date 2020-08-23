
using System;

namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Represents a wavelet defined by an bi-orthogonal base.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Wavelets.OrthogonalWavelet" />
    public class BiorthogonalWavelet : OrthogonalWavelet
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="BiorthogonalWavelet"/> class.
        /// </summary>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The family name.</param>
        /// <param name="lowRec">The low pass reconstruction filter.</param>
        /// <param name="highRec">The high pass reconstruction filter.</param>
        /// <param name="lowDec">The low pass decomposition filter.</param>
        /// <param name="highDec">The high pass decomposition filter.</param>
        /// <param name="vanishingMoments">The vanishing moments.</param>
        /// <param name="freq">The central frequency.</param>
        /// <exception cref="ArgumentNullException">Thrown when any of the filters is <c>null</c>.</exception>
        /// <exception cref="ArgumentException">Thrown when the filter lengths do not match.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="vanishingMoments"/> is negative.</exception>
        public BiorthogonalWavelet(string shortName, string familyName, double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, int vanishingMoments, double freq = 0)
             : base(shortName, familyName, lowRec, highRec, lowDec, highDec, vanishingMoments, freq)
        {
            Other = new BiorthogonalWavelet(this, shortName, familyName, lowDec.Reversed(), highDec.Reversed(), lowRec.Reversed(), highRec.Reversed(), vanishingMoments, freq);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BiorthogonalWavelet"/> class.
        /// </summary>
        /// <param name="reverse">The other wavelet.</param>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The family name.</param>
        /// <param name="lowRec">The low pass reconstruction filter.</param>
        /// <param name="highRec">The high pass reconstruction filter.</param>
        /// <param name="lowDec">The low pass decomposition filter.</param>
        /// <param name="highDec">The high pass decomposition filter.</param>
        /// <param name="vanishingMoments">The vanishing moments.</param>
        /// <param name="freq">The central frequency.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="reverse"/> or any of the filters is <c>null</c>.</exception>
        /// <exception cref="ArgumentException">Thrown when the filter lengths do not match.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="vanishingMoments"/> is negative.</exception>
        protected BiorthogonalWavelet(BiorthogonalWavelet reverse, string shortName, string familyName, double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, int vanishingMoments, double freq = 0)
             : base(shortName, familyName, lowRec, highRec, lowDec, highDec, vanishingMoments, freq)
        {
            Other = reverse ?? throw new ArgumentNullException(nameof(reverse));
        }

        /// <summary>
        /// Gets the other wavelet function.
        /// </summary>
        public BiorthogonalWavelet Other { get; private set; }

        /// <summary>
        /// Gets the size that should have the signal returned by <see cref="M:Neuronic.TimeFrequency.Wavelets.OrthogonalWavelet.EvaluateDomain" />.
        /// </summary>
        protected override int DesiredOutputSize => base.DesiredOutputSize - 1;

        /// <summary>
        /// Convolves the specified coefficient vector with the filters to the specified level of precision.
        /// </summary>
        /// <param name="coeffs">The coefficients to convolve.</param>
        /// <param name="level">The level.</param>
        /// <returns>
        /// The convolved coefficients.
        /// </returns>
        protected override double[] Upcoef(double[] coeffs, int level)
        {
            if (VanishingMoments % 4 != 1)
                coeffs[0] *= -1;
            return base.Upcoef(coeffs, level);
        }
    }
}
