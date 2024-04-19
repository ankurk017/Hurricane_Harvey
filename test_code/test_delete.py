def calculate_error(out, harvey):
 ws_error = []
 mslp_error = []
 track_error = []
 for out1a in out:
  model_ws = np.array([harvey['vmax'][np.where(harvey['date'] == pd.to_datetime(str(val.values)).to_pydatetime())[0]] for val in out1a['mslp']["Time"]], dtype='object')
  ws_error.append((model_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]-(out1a['ws'][np.array([ws.shape[0] != 0 for ws in model_ws])]).values*1.95).mean())
  
  model_mslp = np.array([harvey['mslp'][np.where(harvey['date'] == pd.to_datetime(str(val.values)).to_pydatetime())[0]] for val in out1a['mslp']["Time"]], dtype='object')
  mslp_error.append((model_mslp[np.array([ws.shape[0] != 0 for ws in model_ws])]-(out1a['mslp'][np.array([ws.shape[0] != 0 for ws in model_ws])]).values).mean())
  
  model_lon = np.array([harvey['lon'][np.where(harvey['date'] == pd.to_datetime(str(val.values)).to_pydatetime())[0]] for val in out1a['mslp']["Time"]], dtype='object')
  model2_lon = (model_lon[np.array([ws.shape[0] != 0 for ws in model_ws])]-(np.array(out1a['track_lon'])[np.array([ws.shape[0] != 0 for ws in model_ws])]))
  model_lat = np.array([harvey['lat'][np.where(harvey['date'] == pd.to_datetime(str(val.values)).to_pydatetime())[0]] for val in out1a['mslp']["Time"]], dtype='object')
  model2_lat = (model_lat[np.array([ws.shape[0] != 0 for ws in model_ws])]-(np.array(out1a['track_lat'])[np.array([ws.shape[0] != 0 for ws in model_ws])]))
  
  track_error.append(np.sqrt(np.concatenate(model2_lon**2 + model2_lat**2)).mean()*111.11)

 return {'ws':ws_error, 'mslp':mslp_error, 'track':track_error}



