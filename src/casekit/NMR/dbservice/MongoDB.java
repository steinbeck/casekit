/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package casekit.NMR.dbservice;

import com.mongodb.MongoClient;
import com.mongodb.MongoClientOptions;
import com.mongodb.MongoCredential;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class MongoDB {


    public static MongoClient login(final String mongoUser, final String mongoPassword, final String mongoAuthDB) throws CDKException {
        MongoClient mongo;
        try {
            // Creating a Mongo client   
            mongo = new MongoClient(
                    new ServerAddress("127.0.0.1", 27017),
                    MongoCredential.createCredential(
                            mongoUser,
                            mongoAuthDB,
                            mongoPassword.toCharArray()),
                    MongoClientOptions.builder().build());
            System.out.println("Login to MongoDB was successfull");
            // Accessing the database             
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(Thread.currentThread().getStackTrace()[1].getMethodName() + ": could not connect to MongoDB!");

            return null;
        }

        return mongo;
    }

    public static MongoDatabase getDatabase(final MongoClient mongo, final String mongoDBName){
        return mongo.getDatabase(mongoDBName);
    }
    
    public static MongoCollection<Document> getCollection(final MongoClient mongo, final String mongoDBName, final String mongoDBCollection) {
        final MongoDatabase database = MongoDB.getDatabase(mongo, mongoDBName);
        if (database == null) {
            return null;
        }
        System.out.println("Access to database \"" + mongoDBName + "\" was successfull");
        // Retrieving a collection
        final MongoCollection<Document> collection = database.getCollection(mongoDBCollection);
        System.out.println("Retrieval of collection \"" + mongoDBCollection + "\" was successfull -> size: " + collection.countDocuments());

        return collection;
    }

    public static void logout(final MongoClient mongo) {
        mongo.close();
    }
}
