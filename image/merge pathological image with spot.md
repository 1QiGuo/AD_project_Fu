# function

```{python}
def merge_transformed_RGB_to_spot_and_save(spot_row_in_fullres,
                                          spot_col_in_fullres,
                                          #max_row, max_col,
                                          imagepath,
                                          #X_transformed,
                                          modulescore,
                                          modulescore_colname,
                                          result_path,
                                          sample_name,
                                          #img_type,
                                          #img_folder,
                                          plot_spot_radius
                                          ):
    #using opencv
    ret, cv_object=cv2.imreadmulti(imagepath, [], cv2.IMREAD_ANYCOLOR)
    cv_object_9 = cv_object[8]
    #it is a gray-scale image
    #change gray-scale image into rgb image
    cv_object_9_t1=np.array([cv_object_9,cv_object_9,cv_object_9]).transpose([1,2,0])
    cv_object_9_t2=np.ascontiguousarray(cv_object_9_t1, dtype=np.uint8)
    overlay = cv_object_9_t2.copy()
    #img = np.ones(shape=(max_row + 1, max_col + 1, 3), dtype=np.uint8) * 255
    
    #let modulescore start from 0
    modulescore_scale=modulescore.copy(deep=True)
    min_score=modulescore_scale[modulescore_colname].min()
    modulescore_scale[modulescore_colname]=modulescore_scale[modulescore_colname]+abs(min_score)
    max_score=modulescore_scale[modulescore_colname].max()
    #scale to 0-255
    modulescore_scale[modulescore_colname]=modulescore_scale[modulescore_colname]*255/max_score
    for index in range(len(modulescore_scale[modulescore_colname])):
        cv2.circle(cv_object_9_t2, 
                   (spot_col_in_fullres.iloc[index], spot_row_in_fullres.iloc[index]),
                      plot_spot_radius,
                      color=(int(modulescore_scale[modulescore_colname].iloc[index]), 0, 256),
                      thickness=-1)
    #add transparancy
    alpha = 0.3
    cv_object_f = cv2.addWeighted(cv_object_9_t2, alpha, overlay, 1 - alpha, 0)
    #export
    cv2.imwrite(result_path+sample_name + 'merged.tiff', cv2.cvtColor(cv_object_f, cv2.COLOR_RGB2BGR))
    #add figure legend
    fig, ax = plt.subplots()
    cmap = clr.LinearSegmentedColormap.from_list('custom blue', ['#0000ff','#ff00ff'], N=255)
    img = ax.imshow(cv_object_f, 
                vmin=modulescore[modulescore_colname].min(), 
                vmax=modulescore[modulescore_colname].max(), 
                cmap=cmap)
    cbar = fig.colorbar(img)
    cbar.ax.set_ylabel('Label', rotation=90)
    plt.savefig(result_path+sample_name + 'merged_colorbar.pdf')
    return modulescore_scale
```

# output

## 2-3
```{python }
result_path="/bmbl_data/qiguo/AD/wenjie/result_Apr16"
imagepath="/bmbl_data/qiguo/AD/AD_home/spatial_rawdata/spatial_8/image/stack/B1 2-3 Stack.tif"
tissue_positioin_23 = pd.read_csv('/bmbl_data/qiguo/AD/AD_home/spatial_rawdata/spatial_8/image/2_3/2_3_tissue_position_filtered.csv', 
                                  delimiter = ',',
                                  index_col=0)
modulescore_23 = pd.read_csv('/bmbl_data/qiguo/AD/AD_home/spatial_rawdata/spatial_8/image/2_3/2_3_modulescore.csv', 
                             delimiter = ',',
                             index_col=0)
diameter2_3=math.ceil(232.34944000000002*0.5)                         

#plot 2-3 spot plot using tissue file, json (diameter), and module score
#2_3 geneset1
module_scale=merge_transformed_RGB_to_spot_and_save(tissue_positioin_23["pxl_row_in_fullres"],
                                       tissue_positioin_23["pxl_col_in_fullres"],
                                       imagepath,
                                       modulescore_23,
                                       "Geneset11",
                                       result_path,
                                       "2_3_geneset1",
                                       diameter2_3#should be integer instead of float
                                         )
#2_3 geneset2
module_scale=merge_transformed_RGB_to_spot_and_save(tissue_positioin_23["pxl_row_in_fullres"],
                                       tissue_positioin_23["pxl_col_in_fullres"],
                                       imagepath,
                                       modulescore_23,
                                       "Geneset21",
                                       result_path,
                                       "2_3_geneset2",
                                       diameter2_3#should be integer instead of float
                                         )
#2_3 geneset1                                        
module_scale=merge_transformed_RGB_to_spot_and_save(tissue_positioin_23["pxl_row_in_fullres"],
                                       tissue_positioin_23["pxl_col_in_fullres"],
                                       imagepath,
                                       modulescore_23,
                                       "Geneset31",
                                       result_path,
                                       "2_3_geneset3",
                                       diameter2_3#should be integer instead of float
                                         )
```

## 1-1
