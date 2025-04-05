#pragma once

#include "CoreMinimal.h"
#include "Engine/DataAsset.h"
#include "PatrolRoute.generated.h"

UCLASS(BlueprintType)
class UPatrolRoute : public UDataAsset
{
	GENERATED_BODY()
    
public:
	UPROPERTY(EditAnywhere, BlueprintReadOnly)
	TArray<FVector> PatrolPoints; 
};
