var year = 2004;
var start_time = year+'-01-01';
var end_time = (year+1)+'-01-01';
var pre_time = (year-1)+'-07-02';
var next_time = (year+1)+'-07-03';
var MOD11A1_2year = ee.ImageCollection('MODIS/061/MOD11A1')
                  .filter(ee.Filter.date(pre_time, next_time))
                  .select(['LST_Day_1km', 'LST_Night_1km']);
                  print(MOD11A1_2year)
var MYD11A1_2year = ee.ImageCollection('MODIS/061/MYD11A1')
                  .filter(ee.Filter.date(pre_time, next_time))
                  .select(['LST_Day_1km', 'LST_Night_1km']);
                  // print(MYD11A1)
                  
var daliySATCollection_2year = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
                      .filter(ee.Filter.date(pre_time, next_time))
                      .select('temperature_2m');
                      
var solution = MOD11A1_2year.first().projection().nominalScale().getInfo();
var matchedImages_2year = MOD11A1_2year.combine(MYD11A1_2year);

var matchedImages_preyear = matchedImages_2year.filter(ee.Filter.date(pre_time, start_time));
var matchedImages_curyear = matchedImages_2year.filter(ee.Filter.date(start_time, end_time));
var matchedImages_nextyear = matchedImages_2year.filter(ee.Filter.date(end_time, next_time));

function getDailyLST(image){
  var terraDayLST=image.select('LST_Day_1km').multiply(0.02);
  var terraNightLST=image.select('LST_Night_1km').multiply(0.02);
  var aquaDayLST = image.select('LST_Day_1km_1').multiply(0.02);
  var aquaNightLST = image.select('LST_Night_1km_1').multiply(0.02);
  terraDayLST = terraDayLST.updateMask(terraDayLST.gt(150).and(terraDayLST.lt(400)));
  terraNightLST = terraNightLST.updateMask(terraNightLST.gt(150).and(terraNightLST.lt(400)));
  aquaDayLST = aquaDayLST.updateMask(aquaDayLST.gt(150).and(aquaDayLST.lt(400)));
  aquaNightLST = aquaNightLST.updateMask(aquaNightLST.gt(150).and(aquaNightLST.lt(400)));

    var TD_TN_AD_AN = terraDayLST.multiply(0.1807).add(terraNightLST.multiply(0.3210))
      .add(aquaDayLST.multiply(0.1907)).add(aquaNightLST.multiply(0.3241)).add(-4.75);
    
    TD_TN_AD_AN = TD_TN_AD_AN.updateMask(terraDayLST.mask()
    .and(terraNightLST.mask())
    .and(aquaDayLST.mask())
    .and(aquaNightLST.mask()));
    
    var TN_AD_AN= terraNightLST.multiply(0.3243).add(aquaDayLST.multiply(0.3582))
    .add(aquaNightLST.multiply(0.3318)).add(-4.31);
    TN_AD_AN = TN_AD_AN.updateMask(terraDayLST.mask().not()
    .and(terraNightLST.mask())
    .and(aquaDayLST.mask())
    .and(aquaNightLST.mask()));

    var TD_TN_AN = terraDayLST.multiply(0.3665).add(terraNightLST.multiply(0.3354))
      .add(aquaNightLST.multiply(0.3216)).add(-6.26);
    TD_TN_AN = TD_TN_AN.updateMask(terraDayLST.mask()
    .and(terraNightLST.mask())
    .and(aquaDayLST.mask().not())
    .and(aquaNightLST.mask()));

    var TD_AD_AN = terraDayLST.multiply(0.1942)
      .add(aquaDayLST.multiply(0.2437)).add(aquaNightLST.multiply(0.5528)).add(2.19);
    TD_AD_AN = TD_AD_AN.updateMask(terraDayLST.mask()
    .and(terraNightLST.mask().not())
    .and(aquaDayLST.mask())
    .and(aquaNightLST.mask()));

    var TD_TN_AD = terraDayLST.multiply(0.2172).add(terraNightLST.multiply(0.5875))
      .add(aquaDayLST.multiply(0.1802)).add(2.88);
    TD_TN_AD = TD_TN_AD.updateMask(terraDayLST.mask()
    .and(terraNightLST.mask())
    .and(aquaDayLST.mask())
    .and(aquaNightLST.mask().not()));

    var TN_AD= terraNightLST.multiply(0.5992)
      .add(aquaDayLST.multiply(0.3821)).add(3.64);
    TN_AD = TN_AD.updateMask(terraDayLST.mask().not()
    .and(terraNightLST.mask())
    .and(aquaDayLST.mask())
    .and(aquaNightLST.mask().not()));
    
    var AD_AN = aquaDayLST.multiply(0.4244)
    .add(aquaNightLST.multiply(0.5637)).add(2.75);
    AD_AN = AD_AN.updateMask(terraDayLST.mask().not()
    .and(terraNightLST.mask().not())
    .and(aquaDayLST.mask())
    .and(aquaNightLST.mask()));

    var TD_AN = terraDayLST.multiply(0.4354)
      .add(aquaNightLST.multiply(0.5630)).add(0.64);
    TD_AN = TD_AN.updateMask(terraDayLST.mask()
    .and(terraNightLST.mask().not())
    .and(aquaDayLST.mask().not())
    .and(aquaNightLST.mask()));
    
    var TD_TN= terraDayLST.multiply(0.3925).add(terraNightLST.multiply(0.5993)).add(1.40);
    TD_TN = TD_TN.updateMask(terraDayLST.mask()
    .and(terraNightLST.mask())
    .and(aquaDayLST.mask().not())
    .and(aquaNightLST.mask().not()));
    
    var mergedImage = TD_TN_AD_AN.addBands([TN_AD_AN, TD_TN_AN, TD_AD_AN, TD_TN_AD
    , TN_AD, AD_AN, TD_AN, TD_TN])
    .rename(['TD_TN_AD_AN', 'TN_AD_AN', 'TD_TN_AN', 'TD_AD_AN', 'TD_TN_AD'
    , 'TN_AD', 'AD_AN', 'TD_AN', 'TD_TN']);
    
    
    mergedImage = mergedImage.reduce(ee.Reducer.mean()).rename('LST');
    mergedImage = mergedImage.set('system:time_start',ee.Date(image.get('system:time_start')));
    var doy = ee.Date(mergedImage.get('system:time_start')).getRelative('day', 'year');
    mergedImage = mergedImage.addBands(ee.Image.constant(1))
    return mergedImage;
}

var pre_year = year-1;
var sub_pre_doy;
if ((pre_year % 4 === 0 && pre_year % 100 !== 0) || (pre_year % 400 === 0)) {
  sub_pre_doy = 366;
} else {
  sub_pre_doy = 365;
}

var sub_cur_doy;
if ((year % 4 === 0 && year % 100 !== 0) || (year % 400 === 0)) {
  sub_cur_doy = 366;
} else {
  sub_cur_doy = 365;
}

var dailyLST_preyear = matchedImages_preyear.map(function(image){
  var doy = ee.Date(image.get('system:time_start')).getRelative('day', 'year').subtract(sub_pre_doy);
  var dailyLST = getDailyLST(image);
  dailyLST = dailyLST.addBands(ee.Image(doy).rename('doy'));
  return dailyLST;
});

var dailyLST_curyear = matchedImages_curyear.map(function(image){
  var doy = ee.Date(image.get('system:time_start')).getRelative('day', 'year');
  var dailyLST = getDailyLST(image);
  dailyLST = dailyLST.addBands(ee.Image(doy).rename('doy'));
  return dailyLST;
});

var dailyLST_nextyear = matchedImages_nextyear.map(function(image){
  var doy = ee.Date(image.get('system:time_start')).getRelative('day', 'year').add(sub_cur_doy);
  var dailyLST = getDailyLST(image);
  dailyLST = dailyLST.addBands(ee.Image(doy).rename('doy'));
  return dailyLST;
});

var lstCombinedIamges_2year = ee.ImageCollection(dailyLST_preyear.toList(dailyLST_preyear.size()).cat(dailyLST_curyear.toList(dailyLST_curyear.size())).cat(dailyLST_nextyear.toList(dailyLST_nextyear.size())));


print(lstCombinedIamges_2year.filter(ee.Filter.date(year-1+'-12-31')).first());
print(lstCombinedIamges_2year.filter(ee.Filter.date(year+'-01-01')).first());
print(lstCombinedIamges_2year.filter(ee.Filter.date(year+'-12-31')).first());
print(lstCombinedIamges_2year.filter(ee.Filter.date(year+1+'-01-01')).first());
var proj = lstCombinedIamges_2year.select('LST').first().projection();
var vaildCount_2year = lstCombinedIamges_2year.select('LST').reduce(ee.Reducer.count()).rename('validCount_2year').reproject(proj);

daliySATCollection_2year = daliySATCollection_2year.map(function(image){
  var daliySAT = image.set('system:time_start',ee.Date(image.get('system:time_start')).millis()).rename('SAT');
  // var doy = ee.Date(daliySAT.get('system:time_start')).getRelative('day', 'year');
  daliySAT = daliySAT.addBands(ee.Image.constant(1));                      
  return daliySAT;
});
var daliySATCollectionList_2year = daliySATCollection_2year.toList(daliySATCollection_2year.size());
daliySATCollection_2year = daliySATCollection_2year.map(function(image){
  var index = daliySATCollectionList_2year.indexOf(image);
  var doy = index.subtract(183);
  image = image.addBands(ee.Image(doy).rename('doy'));
  return image;
});
var harmonicSAT_2year = daliySATCollection_2year.map(function(image){
  var timeRadians = image.select('doy').multiply(2 * Math.PI).divide(365.25);
  var timeRadians2 = image.select('doy').multiply(4 * Math.PI).divide(365.25);
  var harmonicSAT = image
    .addBands(timeRadians.cos().rename('cos1'))
    .addBands(timeRadians.sin().rename('sin1'))
    .addBands(timeRadians2.cos().rename('cos2'))
    .addBands(timeRadians2.sin().rename('sin2'));
  return harmonicSAT;
});

var harmonicIndependentsSAT_2year = ee.List(['constant', 'cos1', 'sin1', 'cos2', 'sin2']);

var harmonicTrendSAT_2year = harmonicSAT_2year
  .select(harmonicIndependentsSAT_2year.add('SAT'))
  .reduce(ee.Reducer.linearRegression(harmonicIndependentsSAT_2year.length(), 1));
  
var harmonicTrendCoefficientsSAT_2year = harmonicTrendSAT_2year.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([harmonicIndependentsSAT_2year]);
var fittedHarmonicSAT_2year = harmonicSAT_2year.map(function(image) {
  return image.addBands(
    image.select(harmonicIndependentsSAT_2year)
      .multiply(harmonicTrendCoefficientsSAT_2year)
      .reduce('sum')
      .rename('fitted'));
});


var deltaSAT_2year = fittedHarmonicSAT_2year.map(function(image){
  var differentSAT = image.select('SAT').subtract(image.select('fitted'));
  differentSAT = differentSAT.resample('bicubic').updateMask(differentSAT.mask()).rename('deltaSAT');
  differentSAT = differentSAT.set('system:time_start',ee.Date(image.get('system:time_start')));
  return differentSAT;
  })


  
var deltaSATList_2year = deltaSAT_2year.toList(deltaSAT_2year.size());
deltaSAT_2year = ee.ImageCollection(deltaSATList_2year.map(function(image){
  image = ee.Image(image);
  image = image.set('system:index',ee.Date(image.get('system:time_start')).format("YYYY_MM_dd"));
  return image;
}));

var harmonicLST_2year = lstCombinedIamges_2year.map(function(image){
  var timeRadians = image.select('doy').multiply(2 * Math.PI).divide(365.25);
  var timeRadians2 = image.select('doy').multiply(4 * Math.PI).divide(365.25);
  var harmonicLST = image
    .addBands(timeRadians.cos().rename('cos1'))
    .addBands(timeRadians.sin().rename('sin1'))
    .addBands(timeRadians2.cos().rename('cos2'))
    .addBands(timeRadians2.sin().rename('sin2'));
  return harmonicLST;
});
var harmonicLSTWithSAT_2year = harmonicLST_2year.combine(deltaSAT_2year);
print('harmonicLSTWithSAT_2year', harmonicLSTWithSAT_2year)

var harmonicLSTWithSATReprojected_2year = harmonicLSTWithSAT_2year.map(function(image){
  return image.reproject(proj);
});

var harmonicIndependentsLST_2year = ee.List(['constant', 'cos1', 'sin1', 'cos2', 'sin2','deltaSAT']);

var harmonicTrendLST_2year = harmonicLSTWithSATReprojected_2year
  .select(harmonicIndependentsLST_2year.add('LST'))
  .reduce(ee.Reducer.linearRegression(harmonicIndependentsLST_2year.length(), 1));
  
  
var harmonicTrendCoefficientsLST_2year = harmonicTrendLST_2year.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([harmonicIndependentsLST_2year]);
  
harmonicTrendCoefficientsLST_2year = harmonicTrendCoefficientsLST_2year.reproject(proj)

var amplitude_1 = harmonicTrendCoefficientsLST_2year.select('cos1').pow(2).add(harmonicTrendCoefficientsLST_2year.select('sin1').pow(2)).sqrt()
              .add(harmonicTrendCoefficientsLST_2year.select('cos2').pow(2).add(harmonicTrendCoefficientsLST_2year.select('sin2').pow(2)).sqrt());
var harmonicLSTWithSATMasked_2year = harmonicLSTWithSATReprojected_2year.map(function(image){
  image = image.updateMask(vaildCount_2year.gte(12));
  return image;
});
harmonicLSTWithSATMasked_2year = harmonicLSTWithSATMasked_2year.filter(ee.Filter.date(start_time, end_time));
print('harmonicLSTWithSATMasked_2year', harmonicLSTWithSATMasked_2year)

var fittedHarmonicLST_2year = harmonicLSTWithSATMasked_2year.map(function(image) {
  return image.addBands(
    image.select(harmonicIndependentsLST_2year)
      .multiply(harmonicTrendCoefficientsLST_2year)
      .reduce('sum')
      .rename('ATCE_LST')).multiply(50).toInt16().set('system:time_start',ee.Date(image.get('system:time_start')));
});
print('fittedHarmonicLST_2year',fittedHarmonicLST_2year);
var ATCE = fittedHarmonicLST_2year.select('ATCE_LST').filter(ee.Filter.date(start_time, end_time)).toBands();

print(ATCE);

Export.image.toDrive(
  { image: ATCE,
    description:'',
    folder:'',
    region: '',
    scale: solution,  
    crs:'SR-ORG:6974',
    maxPixels:300000000
  }
);
