// Place an absorber layer inside the calorimeter layer
//
Box AbsLayer(AbsLayerX / 2., AbsLayerY / 2., AbsLayerZ / 2.);
Volume AbsLayerVol(
    "AbsLayerVol", AbsLayer,
    description.material(x_abslayer.attr<std::string>(_U(material))));
AbsLayerVol.setVisAttributes(description, x_abslayer.visStr());
PlacedVolume AbsLayerPlaced = CaloLayerVol.placeVolume(
    AbsLayerVol, 1, Position(0., 0., -CaloLayerZ / 2. + AbsLayerZ / 2.));
AbsLayerPlaced.addPhysVolID("abslayer", 1);

// Place an active layer inside the calorimeter layer
//
Box SensLayer(SensLayerX / 2., SensLayerY / 2., SensLayerZ / 2.);
Volume SensLayerVol(
    "SensLayerVol", SensLayer,
    description.material(x_senslayer.attr<std::string>(_U(material))));
SensLayerVol.setVisAttributes(description, x_senslayer.visStr());
PlacedVolume SensLayerPlaced = CaloLayerVol.placeVolume(
    SensLayerVol, 1, Position(0., 0., CaloLayerZ / 2. - SensLayerZ / 2.));
SensLayerPlaced.addPhysVolID("abslayer", 0);
