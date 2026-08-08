/**
 * Pixel-level quality information for daily MODIS observations and Step 1
 * temporal upscaling in Google Earth Engine.
 *
 *
 */


// =============================================================================
// 1. User-defined parameter
// =============================================================================

var TARGET_YEAR = 2005;


// =============================================================================
// 2. Temporal extent and MODIS data
// =============================================================================

var startDate = ee.Date.fromYMD(TARGET_YEAR, 1, 1);
var endDate = startDate.advance(1, 'year');

var expectedDateNames = ee.List.sequence(0,endDate.difference(startDate, 'day').subtract(1))
  .map(function(dayOffset) {
    return startDate
      .advance(ee.Number(dayOffset), 'day')
      .format('YYYY-MM-dd');
  });

var sourceBandNames = [
  'LST_Day_1km',
  'QC_Day',
  'LST_Night_1km',
  'QC_Night'
];

var terraLSTCollection = ee.ImageCollection('MODIS/061/MOD11A1')
  .filterDate(startDate, endDate)
  .select(sourceBandNames, [
    'LST_Day_1km',
    'QC_Day',
    'LST_Night_1km',
    'QC_Night'
  ]);

var aquaLSTCollection = ee.ImageCollection('MODIS/061/MYD11A1')
  .filterDate(startDate, endDate)
  .select(sourceBandNames, [
    'LST_Day_1km_1',
    'QC_Day_1',
    'LST_Night_1km_1',
    'QC_Night_1'
  ]);

var modisProjection = terraLSTCollection
  .first()
  .select('LST_Day_1km')
  .projection();
var modisNominalScale = modisProjection.nominalScale().getInfo();

// Terra and Aqua images are matched by their common daily system:index.
var matchedMODISCollection = terraLSTCollection
  .combine(aquaLSTCollection, false)
  .sort('system:time_start');


// =============================================================================
// 3. Daily UInt16 quality encoding
// =============================================================================

/**
 * Encodes MODIS observation quality and Step 1 temporal-upscaling validity for
 * one pixel-day. TD, TN, AD, and AN denote Terra daytime, Terra nighttime,
 * Aqua daytime, and Aqua nighttime observations, respectively.
 */
function encodeDailyQualityInformation(image) {
  var acquisitionDate = ee.Date(image.get('system:time_start'));
  var dateName = acquisitionDate.format('YYYY-MM-dd');

  var rawLST = image.select(
    [
      'LST_Day_1km',
      'LST_Night_1km',
      'LST_Day_1km_1',
      'LST_Night_1km_1'
    ],
    ['TD', 'TN', 'AD', 'AN']
  );

  var qualityControl = image.select(
    ['QC_Day', 'QC_Night', 'QC_Day_1', 'QC_Night_1'],
    ['TD', 'TN', 'AD', 'AN']
  );

  var rangeValid = rawLST
    .gt(7500)
    .and(rawLST.lt(20000))
    .unmask(0)
    .toUint16();

  // Mandatory QA values 0 and 1 indicate that an LST value was produced.
  var mandatoryQAValid = qualityControl
    .bitwiseAnd(3)
    .lte(1);

  var observationValid = rangeValid
    .and(mandatoryQAValid)
    .unmask(0)
    .toUint16();

  // QC bits 6-7 store the MODIS LST error flags. Adding 1 reserves zero for
  // invalid observations and maps valid error flags to status values 1-4.
  var observationStatus = qualityControl
    .rightShift(6)
    .bitwiseAnd(3)
    .add(1)
    .multiply(observationValid)
    .unmask(0)
    .toUint16();

  var validObservationCount = observationStatus
    .gt(0)
    .reduce(ee.Reducer.sum())
    .toUint16();

  // Encode the availability of TD, TN, AD, and AN.
  var rangeValidityPattern = rangeValid
    .select('TD')
    .bitwiseOr(rangeValid.select('TN').leftShift(1))
    .bitwiseOr(rangeValid.select('AD').leftShift(2))
    .bitwiseOr(rangeValid.select('AN').leftShift(3));

  // Encode Step 1 validity
  var temporalUpscalingValid = rangeValid
    .reduce(ee.Reducer.sum())
    .gte(3)
    .or(rangeValidityPattern.eq(6))
    .or(rangeValidityPattern.eq(12))
    .or(rangeValidityPattern.eq(9))
    .or(rangeValidityPattern.eq(3))
    .unmask(0)
    .toUint16();

  // Pack the six quality fields into a single UInt16 image.
  var packedQuality = observationStatus
    .select('TD')
    .bitwiseOr(observationStatus.select('TN').leftShift(3))
    .bitwiseOr(observationStatus.select('AD').leftShift(6))
    .bitwiseOr(observationStatus.select('AN').leftShift(9))
    .bitwiseOr(validObservationCount.leftShift(12))
    .bitwiseOr(temporalUpscalingValid.leftShift(15))
    .unmask(0)
    .toUint16()
    .rename('quality_status');

  return packedQuality.set({
    'system:index': dateName,
    'system:time_start': acquisitionDate.millis(),
    date: dateName
  });
}

var dailyQualityCollection = matchedMODISCollection
  .map(encodeDailyQualityInformation)
  .sort('system:time_start');


// =============================================================================
// 4. Annual quality-information image
// =============================================================================

var availableDateNames = ee.List(
  dailyQualityCollection.aggregate_array('date')
);
var missingDateNames = expectedDateNames.removeAll(availableDateNames);

var annualQualityInformation = dailyQualityCollection
  .toBands()
  .rename(availableDateNames)
  .toUint16()
  .set({
    year: TARGET_YEAR,
    band_count: availableDateNames.length(),
    data_type: 'UInt16',
    band_name_format: 'YYYY-MM-dd',
    bit_2_0: 'Terra daytime observation flag',
    bit_5_3: 'Terra nighttime observation flag',
    bit_8_6: 'Aqua daytime observation flag',
    bit_11_9: 'Aqua nighttime observation flag',
    bit_14_12: 'Number of QC-valid MODIS observations',
    bit_15: 'Validity of Step 1 temporal upscaling'
  });



// =============================================================================
// 5. Export
// =============================================================================


Export.image.toDrive({
  image: annualQualityInformation,
  description: 'MODIS_observation_temporal_upscaling_QA_' + TARGET_YEAR,
  folder: 'GEE_exports',
  region: exportRegion,
  scale: modisNominalScale,
  crs: 'SR-ORG:6974',
  maxPixels: 1e13
});
