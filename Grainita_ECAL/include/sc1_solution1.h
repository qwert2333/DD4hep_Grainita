// Create the calorimeter layer that will contain both the absorber
// and active elements
//
Box CaloLayer(CaloLayerX / 2., CaloLayerY / 2., CaloLayerZ / 2.);
Volume CaloLayerVol(
    "CaloLayerVol", CaloLayer,
    description.material(x_calolayer.attr<std::string>(_U(material))));
CaloLayerVol.setVisAttributes(description, x_calolayer.visStr());
