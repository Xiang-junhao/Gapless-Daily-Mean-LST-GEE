/**
 * Daily mean land surface temperature reconstruction using the enhanced
 * annual temperature cycle (EATC) model in Google Earth Engine.
 *
 * Workflow:
 * 1. Estimate gapped daily mean LST from two to four valid Terra and Aqua
 *    MODIS instantaneous observations.
 * 2. Derive the daily ERA5-Land air-temperature residual from ATC model.
 * 3. Fit the EATC model using gapped daily mean LST and the air-temperature residual.
 * 4. Generate gapless daily LST for the target year.
 *
 */


// =============================================================================
// 1. User-defined parameters
// =============================================================================

var TARGET_YEAR = 2004;


// =============================================================================
// 2. Temporal extent and input data
// =============================================================================

// The fitting window extends using a two-year moving window (current year ± six months)
var targetStartDate = ee.Date.fromYMD(TARGET_YEAR, 1, 1);
var targetEndDate = targetStartDate.advance(1, 'year');
var fittingStartDate = ee.Date.fromYMD(TARGET_YEAR - 1, 7, 2);
var fittingEndDate = ee.Date.fromYMD(TARGET_YEAR + 1, 7, 3);

var terraLSTCollection = ee.ImageCollection('MODIS/061/MOD11A1')
  .filterDate(fittingStartDate, fittingEndDate)
  .select(['LST_Day_1km', 'LST_Night_1km']);

var aquaLSTCollection = ee.ImageCollection('MODIS/061/MYD11A1')
  .filterDate(fittingStartDate, fittingEndDate)
  .select(['LST_Day_1km', 'LST_Night_1km']);

var dailyAirTemperatureCollection = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
  .filterDate(fittingStartDate, fittingEndDate)
  .select('temperature_2m');

var modisProjection = terraLSTCollection
  .first()
  .select('LST_Day_1km')
  .projection();
var modisNominalScale = modisProjection.nominalScale().getInfo();

// Terra and Aqua images are matched by their common daily system:index.
var matchedMODISCollection = terraLSTCollection.combine(aquaLSTCollection);


// =============================================================================
// 3. Daily mean LST estimation
// =============================================================================

/**
 * Estimates daily mean LST using the regression coefficients of Xing et al.
 * (2021). TD, TN, AD, and AN denote Terra daytime, Terra nighttime, Aqua
 * daytime, and Aqua nighttime observations, respectively.
 */
function estimateDailyMeanLST(image) {
  // Convert MODIS values to kelvin and retain values within 150-400 K.
  var terraDayLST = image.select('LST_Day_1km').multiply(0.02).rename('TD');
  var terraNightLST = image.select('LST_Night_1km').multiply(0.02).rename('TN');
  var aquaDayLST = image.select('LST_Day_1km_1').multiply(0.02).rename('AD');
  var aquaNightLST = image.select('LST_Night_1km_1').multiply(0.02).rename('AN');

  terraDayLST = terraDayLST.updateMask(
    terraDayLST.gt(150).and(terraDayLST.lt(400))
  );
  terraNightLST = terraNightLST.updateMask(
    terraNightLST.gt(150).and(terraNightLST.lt(400))
  );
  aquaDayLST = aquaDayLST.updateMask(
    aquaDayLST.gt(150).and(aquaDayLST.lt(400))
  );
  aquaNightLST = aquaNightLST.updateMask(
    aquaNightLST.gt(150).and(aquaNightLST.lt(400))
  );

  var tdValid = terraDayLST.mask();
  var tnValid = terraNightLST.mask();
  var adValid = aquaDayLST.mask();
  var anValid = aquaNightLST.mask();

  // Estimate based on all four observations.
  var estimateTDTNADAN = terraDayLST
    .multiply(0.1807)
    .add(terraNightLST.multiply(0.3210))
    .add(aquaDayLST.multiply(0.1907))
    .add(aquaNightLST.multiply(0.3241))
    .add(-4.75)
    .updateMask(tdValid.and(tnValid).and(adValid).and(anValid));

  // Estimates based on three observations.
  var estimateTNADAN = terraNightLST
    .multiply(0.3243)
    .add(aquaDayLST.multiply(0.3582))
    .add(aquaNightLST.multiply(0.3318))
    .add(-4.31)
    .updateMask(tdValid.not().and(tnValid).and(adValid).and(anValid));

  var estimateTDTNAN = terraDayLST
    .multiply(0.3665)
    .add(terraNightLST.multiply(0.3354))
    .add(aquaNightLST.multiply(0.3216))
    .add(-6.26)
    .updateMask(tdValid.and(tnValid).and(adValid.not()).and(anValid));

  var estimateTDADAN = terraDayLST
    .multiply(0.1942)
    .add(aquaDayLST.multiply(0.2437))
    .add(aquaNightLST.multiply(0.5528))
    .add(2.19)
    .updateMask(tdValid.and(tnValid.not()).and(adValid).and(anValid));

  var estimateTDTNAD = terraDayLST
    .multiply(0.2172)
    .add(terraNightLST.multiply(0.5875))
    .add(aquaDayLST.multiply(0.1802))
    .add(2.88)
    .updateMask(tdValid.and(tnValid).and(adValid).and(anValid.not()));

  // Estimates based on two daytime-nighttime observations.
  var estimateTNAD = terraNightLST
    .multiply(0.5992)
    .add(aquaDayLST.multiply(0.3821))
    .add(3.64)
    .updateMask(tdValid.not().and(tnValid).and(adValid).and(anValid.not()));

  var estimateADAN = aquaDayLST
    .multiply(0.4244)
    .add(aquaNightLST.multiply(0.5637))
    .add(2.75)
    .updateMask(tdValid.not().and(tnValid.not()).and(adValid).and(anValid));

  var estimateTDAN = terraDayLST
    .multiply(0.4354)
    .add(aquaNightLST.multiply(0.5630))
    .add(0.64)
    .updateMask(tdValid.and(tnValid.not()).and(adValid.not()).and(anValid));

  var estimateTDTN = terraDayLST
    .multiply(0.3925)
    .add(terraNightLST.multiply(0.5993))
    .add(1.40)
    .updateMask(tdValid.and(tnValid).and(adValid.not()).and(anValid.not()));

  var dailyMeanLST = estimateTDTNADAN
    .addBands([
      estimateTNADAN,
      estimateTDTNAN,
      estimateTDADAN,
      estimateTDTNAD,
      estimateTNAD,
      estimateADAN,
      estimateTDAN,
      estimateTDTN
    ])
    .reduce(ee.Reducer.mean())
    .rename('LST');

  var acquisitionDate = ee.Date(image.get('system:time_start'));
  return dailyMeanLST
    .set('system:time_start', acquisitionDate.millis())
    .set('system:index', acquisitionDate.format('YYYY_MM_dd'));
}

// Add Day of year (DOY) band
var dailyMeanLSTCollection = matchedMODISCollection.map(function(image) {
  var dailyMeanLST = estimateDailyMeanLST(image);
  var relativeDay = ee.Date(dailyMeanLST.get('system:time_start'))
    .difference(targetStartDate, 'day');

  return dailyMeanLST.addBands(
    ee.Image.constant(relativeDay).rename('doy').toFloat()
  );
});

var validDailyMeanLSTCount = dailyMeanLSTCollection
  .select('LST')
  .reduce(ee.Reducer.count())
  .rename('valid_count')
  .reproject(modisProjection);


// =============================================================================
// 4. Air-temperature ATC residual
// =============================================================================

// Add two harmonic terms.
function addAnnualHarmonicTerms(image) {
  var firstHarmonicAngle = image
    .select('doy')
    .multiply(2 * Math.PI)
    .divide(365.25);
  var secondHarmonicAngle = image
    .select('doy')
    .multiply(4 * Math.PI)
    .divide(365.25);

  return image
    .addBands(ee.Image.constant(1).rename('constant').toFloat())
    .addBands(firstHarmonicAngle.cos().rename('cos1'))
    .addBands(firstHarmonicAngle.sin().rename('sin1'))
    .addBands(secondHarmonicAngle.cos().rename('cos2'))
    .addBands(secondHarmonicAngle.sin().rename('sin2'));
}

var preparedAirTemperatureCollection = dailyAirTemperatureCollection.map(
  function(image) {
    var acquisitionDate = ee.Date(image.get('system:time_start'));
    var relativeDay = acquisitionDate.difference(targetStartDate, 'day');

    return image
      .select('temperature_2m')
      .rename('SAT')
      .toFloat()
      .addBands(ee.Image.constant(relativeDay).rename('doy').toFloat())
      .set('system:time_start', acquisitionDate.millis())
      .set('system:index', acquisitionDate.format('YYYY_MM_dd'));
  }
);

var harmonicAirTemperatureCollection = preparedAirTemperatureCollection.map(
  addAnnualHarmonicTerms
);

var AIR_TEMPERATURE_PREDICTORS = ee.List([
  'constant',
  'cos1',
  'sin1',
  'cos2',
  'sin2'
]);

// Fit the annual cycle of ERA5-Land air temperature.
var airTemperatureRegression = harmonicAirTemperatureCollection
  .select(AIR_TEMPERATURE_PREDICTORS.add('SAT'))
  .reduce(
    ee.Reducer.linearRegression(AIR_TEMPERATURE_PREDICTORS.length(), 1)
  );

var airTemperatureCoefficients = airTemperatureRegression
  .select('coefficients')
  .arrayProject([0])
  .arrayFlatten([AIR_TEMPERATURE_PREDICTORS]);

// Calculate daily air-temperature residuals
var airTemperatureResidualCollection = harmonicAirTemperatureCollection.map(
  function(image) {
    var fittedAirTemperature = image
      .select(AIR_TEMPERATURE_PREDICTORS)
      .multiply(airTemperatureCoefficients)
      .reduce(ee.Reducer.sum());

    return image
      .select('SAT')
      .subtract(fittedAirTemperature)
      .rename('deltaSAT')
      .resample('bicubic')
      .set('system:time_start', image.get('system:time_start'))
      .set('system:index', image.get('system:index'));
  }
);


// =============================================================================
// 5. EATC fitting and prediction
// =============================================================================

var harmonicDailyMeanLSTCollection = dailyMeanLSTCollection.map(
  addAnnualHarmonicTerms
);

// Match daily mean LST with the ERA5-Land residual using the common date index.
var eatcInputCollection = harmonicDailyMeanLSTCollection
  .combine(airTemperatureResidualCollection)
  .map(function(image) {
    return image.reproject(modisProjection);
  });

var sufficientObservationMask = validDailyMeanLSTCount.gte(12);
var maskedEATCInputCollection = eatcInputCollection.map(function(image) {
  return image.updateMask(sufficientObservationMask);
});

var EATC_PREDICTORS = ee.List([
  'constant',
  'cos1',
  'sin1',
  'cos2',
  'sin2',
  'deltaSAT'
]);

// Estimate pixel-wise EATC coefficients over the complete fitting window.
var eatcRegression = maskedEATCInputCollection
  .select(EATC_PREDICTORS.add('LST'))
  .reduce(ee.Reducer.linearRegression(EATC_PREDICTORS.length(), 1));

var eatcCoefficients = eatcRegression
  .select('coefficients')
  .arrayProject([0])
  .arrayFlatten([EATC_PREDICTORS])
  .reproject(modisProjection);

var targetYearEATCInputCollection = maskedEATCInputCollection.filterDate(
  targetStartDate,
  targetEndDate
);

// Generate and scale the fitted daily EATC LST values.
var fittedDailyEATCCollection = targetYearEATCInputCollection.map(
  function(image) {
    var fittedLST = image
      .select(EATC_PREDICTORS)
      .multiply(eatcCoefficients)
      .reduce(ee.Reducer.sum())
      .rename('EATC_LST')
      .multiply(50)
      .toInt16();

    return fittedLST
      .set('system:time_start', image.get('system:time_start'))
      .set('system:index', image.get('system:index'));
  }
);

var annualEATCLST = fittedDailyEATCCollection
  .sort('system:time_start')
  .toBands()
  .set({
    target_year: TARGET_YEAR,
    units: 'K',
    scale_factor: 0.02,
  });

print('Annual EATC-fitted LST', annualEATCLST);


// =============================================================================
// 6. Export
// =============================================================================

Export.image.toDrive({
  image: annualEATCLST,
  description: 'EATC_daily_mean_LST_' + TARGET_YEAR,
  folder: 'GEE_exports',
  region: exportRegion,
  scale: modisNominalScale,
  crs: 'SR-ORG:6974',
  maxPixels: 1e13
});
