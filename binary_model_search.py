import os
from glob import glob
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import Binarizer
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from mfe_plots import loadProfileData
from data_helpers import getSpeciesProperty
from endosymbionys import isEndosymbiont


def getProfilesForAllGroups(args):
    files = [[x for x in glob(profilesGlob) if os.path.exists(x)] for profilesGlob in args.use_profile_data]

    profiles = []

    for group, files in enumerate(files):
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics, _1, _2) = loadProfileData(files)
        #self.xdata.append( xdata )
        #self.ydata.append( ydata )
        #self.ydata_nativeonly.append( ydata_nativeonly )
        #self.ydata_shuffleonly.append( ydata_shuffledonly )
        #self.label.append( labels )
        #self.groups.append( groups )
        #self.filesUsed.append( filesUsed )
        print("{} -> {} profiles".format(group, len(biasProfiles)))

        profiles.append(biasProfiles)

    return profiles



def getTraits( taxIds, traits=(("gc-content", "float"), ("ENc-prime", "float"), ("optimum-temperature", "float"), ("is-endosymbiont", "int"), ("is-high-temp", "int")) ):
    df = pd.DataFrame(dict([(x, pd.Series(dtype=t)) for x,t in traits]), index=taxIds)
    
    for trait, _ in traits:
        for taxId in taxIds:
            traitVal = None
            
            if trait=="is-endosymbiont":
                traitVal = isEndosymbiont( taxId )

            elif trait=="is-high-temp":
                traitVal = 0
                prop = getSpeciesProperty(taxId, "optimum-temperature")
                
                if not prop[0] is None:
                    tempVal = float(prop[0])
                    print("{} -> {}".format(taxId, tempVal))
                    if tempVal > 75.0:
                        traitVal = 1
                        print(taxId)
                        
            else:
                prop = getSpeciesProperty(taxId, trait)
                
                if not prop[0] is None:
                    traitVal = float(prop[0])
                    
            if not traitVal is None:
                print("{} {} -> {}".format(taxId, trait, traitVal))
                df.loc[taxId, trait] = traitVal
    return df
            

def addProfileAbsMeanValue( traits, profileData, fromPos, toPos ):

    newDf = pd.DataFrame({"profile-0-abs-range-mean":pd.Series(dtype=float)}, index=profileData[0].keys())
    
    for taxId in profileData[0].keys():
        profile = profileData[0][taxId]
        vals = profile.loc[fromPos:toPos]
        #N = vals.size
        meanVal = vals.mean()
        newDf.loc[taxId] = abs(meanVal)

    return  pd.merge( traits, newDf, left_index=True, right_index=True )
        
def addProfileStdDev( traits, profileData, _1, _2 ):

    newDf = pd.DataFrame({"profile-0-sd":pd.Series(dtype=float)}, index=profileData[0].keys())
    
    for taxId in profileData[0].keys():
        profile = profileData[0][taxId]
        newDf.loc[taxId] = profile.std()

    return  pd.merge( traits, newDf, left_index=True, right_index=True )
        
        
    


def binaryModelSearch( args, traits, range2, trait1, trait2 ):

    # Binarize the response trait
    trait1bin = "{}.bin".format(trait1)
    trait1_thres = 0.14 # 0.17 # 0.14
    s = pd.Series( traits.loc[:,trait1] < trait1_thres, name=trait1bin )
    traits = pd.concat( (traits, s), axis=1 )
    for taxId in traits.index.ravel():
        assert( s.loc[taxId] == traits.loc[taxId,trait1bin] )
    trueClassification = 1

    #traits["optimum-temperature"] = map( lambda x:200 - 20*x, traits["profile-0-abs-range-mean"] )
    print(traits)

    # Remove NAs
    #traitsNoNAs = traits[ ~traits.loc[ :, trait2 ].isna() ]
    traitsNoNAs = traits[ ~traits.loc[ :, trait2 ].isnull() ]
    #print(traitsNoNAs.loc[:,trait1bin])
    #print(sum(traitsNoNAs.loc[:,trait1bin]))
    #print( np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1) )
    #print( traitsNoNAs.loc[:,trait1bin] )

    inputData = np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1)
    responseData = traitsNoNAs.loc[:,trait1bin].reshape(-1,1)

    avgprData = ([], [])
    o = []
    for t in np.linspace( range2[0], range2[1], 100):
        pipeline = make_pipeline( Binarizer(threshold=t), LogisticRegression() )

        print(inputData.shape)
        print(responseData.shape)
        pipeline.fit( inputData, responseData )

        #print("////"*10)
        #print(t)
        #print(pipeline)
        #print(pipeline.steps[1])
        #print("{} + {}*x1".format( pipeline.steps[1][1].intercept_, pipeline.steps[1][1].coef_[0] ) )
        #print( pipeline.predict( np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1) ) )
        #print( pipeline.steps[1][1].predict_proba( np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1) ) )
        averagePR = average_precision_score( responseData, pipeline.predict_proba( inputData )[:,trueClassification] )
        AUC = roc_auc_score( responseData, pipeline.predict_proba( inputData )[:,trueClassification] )
        avgprData[0].append(t)
        #avgprData[1].append(averagePR)
        avgprData[1].append(AUC)


        if( len(avgprData[0]) % 20 == 11 ):

            precision, recall, _ = precision_recall_curve( responseData, pipeline.predict_proba( inputData )[:,trueClassification] )
            o.append("T: {:.3}".format(t))
            o.append("avg. pr: {:.4}".format(averagePR) )
            o.append(str(recall))
            o.append(str(precision))
            print( averagePR )
            plt.step( recall, precision )
            print(recall)
            print(precision)
            plt.title("PR ({:.3})".format(t))
            plt.savefig( "test_precision_recall_{:.2}.pdf".format(t) )
            plt.close()

            ty1, ty2, _ = roc_curve( responseData, pipeline.predict_proba( inputData )[:,trueClassification] )
            plt.step( ty1, ty2 )
            plt.title("ROC ({:.3})".format(t))
            plt.savefig( "test_precision_recall_roc_{:.2}.pdf".format(t) )
            plt.close()
            
            
            plt.scatter( inputData, responseData )
            plt.title("data ({:.3})".format(t))
            plt.savefig( "test_precision_recall_{:.2}_data.pdf".format(t) )
            plt.close()

            #plt.scatter( traitsNoNAs.loc[:,trait1bin], pipeline.steps[1][1].predict_proba( np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1) )[:,0] )
            xx = np.linspace(-1, 4, 500).reshape(-1,1)
            plt.scatter( xx, pipeline.steps[1][1].predict_proba( xx )[:,trueClassification] )
            plt.title("probit ({:.3}c)".format(t))
            plt.savefig( "test_precision_recall_{:.2}_probit.pdf".format(t) )
            plt.close()

            xx = np.linspace(range2[0], range2[1], 500).reshape(-1,1)
            plt.scatter( xx, pipeline.predict_proba( xx )[:,trueClassification] )
            plt.title("pipeline.probit ({:.3}c)".format(t))
            plt.savefig( "test_precision_recall_{:.2}_pipeline_probit.pdf".format(t) )
            plt.close()
            
            
            
        plt.scatter( avgprData[0], avgprData[1] )
        plt.title("Avg. PR by threshold")
        plt.savefig( "test_precision_recall_avgpr_by_thres.pdf" )
        plt.close()

            
    for i in o:
        print(i)
    
    return 
    # Create pipeline
    pipeline = make_pipeline( Binarizer(), LogisticRegression() )
    # Pipeline(memory=None,
    #    steps=[('binarizer', Binarizer(copy=True, threshold=0.0)), ('logisticregression', LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
    #         intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
    #         penalty='l2', random_state=None, solver='liblinear', tol=0.0001,
    #         verbose=0, warm_start=False))])


    # Initialize parameters space to be searched
    #r1 = (traits.loc[:,trait1].min(), traits.loc[:,trait1].max())
    #span_r1 = r1[1]-r1[0]
    r2 = (traits.loc[:,trait2].min(), traits.loc[:,trait2].max())
    span_r2 = r2[1]-r2[0]
    Nsteps = 93
    margin = 1/(Nsteps*2.0)
    
    #params_grid = [{trait1: np.linspace( r1[0]+span_r1*margin, r1[1]-span_r1*margin, Nsteps ),
    #                trait2: np.linspace( r2[0]+span_r2*margin, r2[1]-span_r2*margin, Nsteps ) }]
    params_grid = [{"binarizer__threshold": np.linspace( r2[0]+span_r2*margin, r2[1]-span_r2*margin, Nsteps ) }]
    print(params_grid)
  
    # Initialize grid search
    gs = GridSearchCV( pipeline, params_grid, cv=5 )
    gs.fit( np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1), traitsNoNAs.loc[:,trait1bin] )
    print(gs.best_params_)


def SVMtest( args, traits, trait1, trait2 ):

    # Binarize the response trait
    trait1bin = "{}.bin".format(trait1)
    trait1_thres = 0.14 # 0.17 # 0.14
    s = pd.Series( traits.loc[:,trait1] < trait1_thres, name=trait1bin )
    traits = pd.concat( (traits, s), axis=1 )
    for taxId in traits.index.ravel():
        assert( s.loc[taxId] == traits.loc[taxId,trait1bin] )

    # Remove NAs
    #traitsNoNAs = traits[ ~traits.loc[ :, trait2 ].isna() ]
    #traitsNoNAs = traits[ ~traits.loc[ :, trait2 ].isnull() ]
    #print(traitsNoNAs.loc[:,trait1bin])
    #print(sum(traitsNoNAs.loc[:,trait1bin]))
    #print( np.array(traitsNoNAs.loc[:,trait2]).reshape(-1,1) )
    #print( traitsNoNAs.loc[:,trait1bin] )

    trait2n = len(trait2)
    inputData = traits.loc[:,trait2]
    print(inputData)
    responseData = traits.loc[:,trait1bin].ravel()
    print(inputData.shape)
    print(responseData.shape)
    
    sv = SVC(tol=1e-6, kernel="linear", probability=True )
    sv.fit( inputData, responseData )

    print( dict( zip( trait2, sv.coef_.ravel() ) ) )
    print( sv.score( inputData, responseData ) )

    #AvPrecision = average_precision_score( responseData, sv.decision_function( inputData ) )
    AvPrecision = average_precision_score( responseData, sv.predict_proba( inputData )[:,1] )
    print("Av. precision: {}".format( AvPrecision ) )
    #precision, recall, _ = precision_recall_curve( responseData, sv.decision_function( inputData ) )
    precision, recall, _ = precision_recall_curve( responseData, sv.predict_proba( inputData )[:,1] )
    plt.step( recall, precision, color="blue", where="post" )
    plt.xlabel( "Recall" )
    plt.ylabel( "Precision" )
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.title( "AUPR={:.3}".format(AvPrecision) )
    plt.savefig( "test_precision_recall_svm_aupr.pdf" )
    plt.close()

    AUC = roc_auc_score( responseData, sv.predict_proba( inputData )[:,1] )
    ty1, ty2, _ = roc_curve( responseData, sv.predict_proba( inputData )[:,1] )
    plt.step( ty1, ty2 )
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.title("ROC (AUC={:.3})".format(AUC))
    plt.savefig( "test_precision_recall_svm_roc.pdf" )
    plt.close()

    #a = np.zeros( (len(range(20,80,5)), len(range(40,60,2))) )
    a = np.zeros( (len(range(10,105,5)), len(range(20,90,5))) )
    for m, gc in enumerate(range(10,105,5)):
        for n, enc in enumerate(range(20,90,5)):
            #trait2=("gc-content", "ENc-prime", "is-high-temp", "is-endosymbiont") )
            a[m,n] = sv.predict_proba( np.array( (float(gc), float(enc), 0, 0) ).reshape(1,-1) )[:,1]
    print(a)
    fig, ax = plt.subplots()
    ax.imshow(a)
    plt.savefig("test_precision_recall_svm_space.pdf" )
    plt.close()
            
    fig, ax = plt.subplots()
    ax.scatter( traits.loc[:,"gc-content"], traits.loc[:,"ENc-prime"], c=traits.loc[:,trait1bin],                                           s=5.0, zorder=1, alpha=0.5, edgecolors='none' )
    ax.scatter( traits.loc[:,"gc-content"], traits.loc[:,"ENc-prime"], c=((traits.loc[:,"is-high-temp"]>0.5) | (traits.loc[:,"is-endosymbiont"]>0.5)) , s=2.0, zorder=2, alpha=1.0, edgecolors='none' )
    plt.title( "Weak dLFE" )
    plt.savefig("test_precision_recall_svm_scatter.pdf" )
    plt.close()


    fig, ax = plt.subplots()
    # for taxId in traits.index.ravel():
    #     jin = inputData.loc[:,trait2].loc[taxId]
    #     jout = sv.predict_proba( jin )[:,1]
    #     print("{} -> {}".format(jin, jout))
    #     ax.scatter( traits.loc[taxId,"gc-content"], traits.loc[taxId,"ENc-prime"], c=jout, cmap=plt.cm.plasma_r, s=5.0, zorder=2, alpha=1.0, edgecolors='none' )
    pf = ax.scatter( traits.loc[taxId,"gc-content"], traits.loc[taxId,"ENc-prime"], c=sv.predict_proba( inputData.loc[:,trait2] )[:,1], cmap=plt.cm.plasma_r, s=5.0, zorder=2, alpha=1.0, edgecolors='none' )
    plt.title( "SVM decision" )
    plt.colorbar(pf)
    plt.savefig("test_precision_recall_svm_decision.pdf" )
    plt.close()


def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

def standalone():
    import argparse
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--use-profile-data", type=parseList(str), required=True)
    args = argsParser.parse_args()

    profiles = getProfilesForAllGroups(args)

    traits = addProfileStdDev(
        getTraits( profiles[0].keys() ),
        profiles,
        150, 290 )

    print(traits)

    #return binaryModelSearch( args, traits, (40,110), trait1="profile-0-abs-range-mean", trait2="optimum-temperature" ) # 0.14: 61
    #return binaryModelSearch( args, traits, (20,60), trait1="profile-0-abs-range-mean",  trait2="gc-content" ) # 0.14: 35-54
    #return binaryModelSearch( args, traits, (40,60), trait1="profile-0-abs-range-mean", trait2="ENc-prime" ) # 0.14: 55
    #return SVMtest( args, traits, trait1="profile-0-abs-range-mean", trait2=("gc-content", "ENc-prime", "optimum-temperature", "is-endosymbiont") )
    return SVMtest( args, traits, trait1="profile-0-abs-range-mean", trait2=("gc-content", "ENc-prime", "is-high-temp", "is-endosymbiont") )
        
    return 0


if __name__=="__main__":
    import sys
    sys.exit(standalone())
        

