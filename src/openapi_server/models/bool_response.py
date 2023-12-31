# Coding: utf-8

from __future__ import annotations
from datetime import date, datetime  # Noqa: F401

import re  # Noqa: F401
from typing import Any, Dict, List, Optional  # Noqa: F401

from pydantic import AnyUrl, BaseModel, EmailStr, Field, validator  # Noqa: F401


class BoolResponse(BaseModel):
    """NOTE: This class is auto generated by OpenAPI Generator (https://openapi-generator.tech).

    Do not edit the class manually.

    BoolResponse - a model defined in OpenAPI

        Value: The value of this BoolResponse [Optional].
        ClientTransactionID: The client_transaction_id of this BoolResponse [Optional].
        ServerTransactionID: The server_transaction_id of this BoolResponse [Optional].
        ErrorNumber: The error_number of this BoolResponse [Optional].
        ErrorMessage: The error_message of this BoolResponse [Optional].
    """

    Value: Optional[bool] = Field(alias="Value", default=None)
    ClientTransactionID: Optional[int] = Field(alias="ClientTransactionID", default=None)
    ServerTransactionID: Optional[int] = Field(alias="ServerTransactionID", default=None)
    ErrorNumber: Optional[int] = Field(alias="ErrorNumber", default=None)
    ErrorMessage: Optional[str] = Field(alias="ErrorMessage", default=None)

    @validator("ClientTransactionID")
    def client_transaction_id_max(cls, value):
        assert value <= 4294967295
        return value

    @validator("ClientTransactionID")
    def client_transaction_id_min(cls, value):
        assert value >= 0
        return value

    @validator("ServerTransactionID")
    def server_transaction_id_max(cls, value):
        assert value <= 4294967295
        return value

    @validator("ServerTransactionID")
    def server_transaction_id_min(cls, value):
        assert value >= 0
        return value

    @validator("ErrorNumber")
    def error_number_max(cls, value):
        assert value <= 2147483647
        return value

    @validator("ErrorNumber")
    def error_number_min(cls, value):
        assert value >= -2147483648
        return value

BoolResponse.update_forward_refs()
