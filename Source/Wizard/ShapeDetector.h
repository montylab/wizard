// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "ShapeDetector.generated.h"

/**
 * 
 */
UCLASS()
class WIZARD_API UShapeDetector : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
	
	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static FString DetectShape(); // 1 circle, 2 square, 3 triangle, 0 - can't predict

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static float getCircle2fResemblance(const TArray<FVector2f>& points);
	
	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static bool IsCircle2f(const TArray<FVector2f>& points, float precision = 0.2f);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static float getCircle3fResemblance(const TArray<FVector3f>& points);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static FPlane ComputeBestFitPlane(const TArray<FVector3f>& points);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static TArray<FVector2f> ProjectPointsToPlane(const TArray<FVector3f>& points, const FPlane& plane);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static bool IsPointNearLine(const FVector2f& p, const FVector2f& p1, const FVector2f& p2, float epsilon);


	// Main API
	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static bool IsCircle3f(const TArray<FVector3f>& points, float precision = 0.2f);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static bool IsTriangle(const TArray<FVector2f>& points, float epsilon, int maxIterations = 10);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static TArray<FVector2f> GetPointsCommonLines(const TArray<FVector2f>& points, float epsilon = 100, int maxIterations = 20);

	UFUNCTION(BlueprintCallable, Category = ShapeDetection)
	static TArray<FVector3f> FilterPointsFromSameOrigin(const TArray<FVector3f>& points, float radius = 3, int count = 3);

};
